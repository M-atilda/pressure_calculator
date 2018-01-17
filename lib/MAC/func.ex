#file   func.ex
#author mi-na
#date   18/01/14
#brief  calculate pressure field
#       right side is derived using central difference
#       in recursion, it applies 2 level multi-grid scheme(v cycle) and SOR method
defmodule MAC.Func do


  def derivePre velocitys_field, pressure, bc_field,
    %{:x_size => x_size,
      :y_size => y_size,
      :z_size => z_size}=information,
    %{:max_ite_times => max_ite_times,
      :error_p => error_p,
      :omega => omega,
      :max_res_ratio => max_res_ratio} do
    IO.puts "[Info] start calc right_side #{inspect DateTime.utc_now}"
    right_side = Enum.map(Enum.to_list(0..(z_size-1)), fn(k) ->
      for j <- 0..(y_size-1) do
        for i <- 0..(x_size-1) do
          calcRSide {i,j,k}, velocitys_field, information
        end
      end end)
    IO.puts "[Info] start derivePreRecurse #{inspect DateTime.utc_now}"
    derivePreRecurse(0, right_side, pressure, bc_field,
      information,
      max_ite_times, error_p, omega,
      max_res_ratio, 1, 1)
  end
  def derivePreRecurse ite_times, right_side, pressure, bc_field,
    %{:dx => dx,
      :dy => dy,
      :dz => dz,
      :x_size => x_size,
      :y_size => y_size,
      :z_size => z_size}=information,
    max_ite_times, error_p, omega,
    max_res_ratio, res_value, res_ratio do
    if ite_times > max_ite_times do
      {:bad, pressure}
    else
      {has_calced, {new_pressure, new_res_field, new_res_value}} = try do
                                                                     {true, derivePreStep(pressure, right_side, bc_field, dx,dy,dz,x_size,y_size,z_size, omega)}
                                                                   rescue
                                                                     err ->
                                                                       IO.puts "[Error] #{inspect err}"
                                                                     {false, {pressure, nil, 0}}
                                                                   end
      if has_calced do
        IO.puts "#{round(:math.log(new_res_value) / :math.log(error_p) * 100)}% #{inspect DateTime.utc_now}"
        if new_res_value < error_p do
          {:ok, new_pressure}
        else
          if (ite_times < 3) || (new_res_value < (max_res_ratio * res_ratio * res_value)) do
            derivePreRecurse(ite_times+1, right_side, new_pressure, bc_field,
              information,
              max_ite_times, error_p, omega,
              max_res_ratio, new_res_value, new_res_value / res_value)
          else
            #NOTE: calculate modification equation in rough grid
            restrict_res_field = restrictField new_res_field, information
            modified_pre_field = calcModField restrict_res_field, information, max_ite_times, error_p
            dp_field = extendField modified_pre_field, information
            modified_new_pressure = Enum.map(Enum.to_list(0..(z_size-1)), fn(k) ->
              for j <- 0..(y_size-1) do
                for i <- 0..(x_size-1) do
                  id(new_pressure, {i,j,k}) + id(dp_field, {i,j,k})
                end end end)
            derivePreRecurse(ite_times+1, right_side, modified_new_pressure, bc_field,
              information,
              max_ite_times, error_p, omega,
              max_res_ratio, new_res_value, new_res_value / res_value)
          end
        end
      else
        {:error, pressure}
      end
    end
  end

  def derivePreStep pressure, right_side, bc_field,
    dx, dy, dz,
    x_size, y_size, z_size,
    omega do
    divide_val = 2*(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz))
    dp_field = Enum.map(Enum.to_list(0..(z_size-1)), fn(k) ->
      for j <- 0..(y_size-1) do
        for i <- 0..(x_size-1) do
          if id(bc_field, {i,j,k}) == nil do
            if 0<i && 0<j && 0<k && i<(x_size-1) && j<(y_size-1) && k<(z_size-1) do
              (((id(pressure, {i+1,j,k}) + id(pressure, {i-1,j,k})) / (2*dx)) + ((id(pressure, {i,j+1,k}) + id(pressure, {i,j-1,k})) / (2*dy)) + ((id(pressure, {i,j,k+1}) + id(pressure, {i,j,k-1})) / (2*dz)) - id(right_side, {i,j,k})) / divide_val - id(pressure, {i,j,k})
            else
              min_i = max 0, i-1
              max_i = min (x_size-1), i+1
              min_j = max 0, j-1
              max_j = min (y_size-1), j+1
              min_k = max 0, k-1
              max_k = min (z_size-1), k+1
              x_width = dx * (max_i - min_i)
              y_width = dy * (max_j - min_j)
              z_width = dz * (max_k - min_k)
              (((id(pressure, {max_i,j,k}) + id(pressure, {min_i,j,k})) / x_width) + ((id(pressure, {i,max_j,k}) + id(pressure, {i,min_j,k})) / y_width) + ((id(pressure, {i,j,max_k}) + id(pressure, {i,j,min_k})) / z_width) - id(right_side, {i,j,k})) / divide_val - id(pressure, {i,j,k})
            end
          else
            0
          end
        end
      end end)
  {Enum.map(Enum.to_list(0..(z_size-1)), fn(k) ->
      for j <- 0..(y_size-1) do
        for i <- 0..(x_size-1) do
          if id(bc_field, {i,j,k}) == nil do
            id(pressure, {i,j,k}) + omega * id(dp_field, {i,j,k})
          else
            id(bc_field, {i,j,k})
          end
        end
      end
    end),
   dp_field,
   List.flatten(dp_field)
   |> Enum.map(&(&1*&1))
   |> :lists.sum}
  end
  
  def calcModField residual,
    %{:dx => small_dx,
      :dy => small_dy,
      :dz => small_dz,
      :x_size => small_x_size,
      :y_size => small_y_size,
      :z_size => small_z_size}, max_ite_times, error_p do
    IO.puts "[Info] start calcModField #{inspect DateTime.utc_now}"    
    dx = 2 * small_dx
    dy = 2 * small_dy
    dz = 2 * small_dz
    x_size = round((small_x_size - 1) / 2)
    y_size = round((small_y_size - 1) / 2)
    z_size = round((small_z_size - 1) / 2)
    divide_val = 2*(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz))
    zero_pressure = Enum.map(Enum.to_list(0..(z_size-1)), fn(_k) ->
      for _j <- 0..(y_size-1) do
        for _i <- 0..(x_size-1) do
          0
        end end end)
    {status, new_pre_field} = try do
                               calcModRecurse(0, zero_pressure, residual,
                                 dx, dy, dz,
                                 x_size, y_size, z_size,
                                 :math.sqrt(max_ite_times), :math.sqrt(error_p), divide_val)
                             rescue
                               error ->
                                 IO.puts "[Error] #{inspect error}"
                               {:error, nil}
                             end
    if status == :error do
      zero_pressure
    else
      new_pre_field
    end
  end
  defp calcModRecurse ite_times, pressure, residual,
    dx, dy, dz,
    x_size, y_size, z_size,
    max_ite_times, error_p, divide_val do
    dp_field = Enum.map(Enum.to_list(0..(z_size-1)), fn(k) ->
      for j <- 0..(y_size-1) do
        for i <- 0..(x_size-1) do
          if 0<i && 0<j && 0<k && i<(x_size-1) && j<(y_size-1) && k<(z_size-1) do
            ((((id(pressure, {i+1,j,k}) + id(pressure, {i-1,j,k})) / (2*dx)) + ((id(pressure, {i,j+1,k}) + id(pressure, {i,j-1,k})) / (2*dy)) + ((id(pressure, {i,j,k+1}) + id(pressure, {i,j,k-1})) / (2*dz))) - id(residual, {i,j,k})) / divide_val
          else
            min_i = max 0, i-1
            max_i = min (x_size-1), i+1
            min_j = max 0, j-1
            max_j = min (y_size-1), j+1
            min_k = max 0, k-1
            max_k = min (z_size-1), k+1
            x_width = dx * (max_i - min_i)
            y_width = dy * (max_j - min_j)
            z_width = dz * (max_k - min_k)
            ((((id(pressure, {max_i,j,k}) + id(pressure, {min_i,j,k})) / x_width) + ((id(pressure, {i,max_j,k}) + id(pressure, {i,min_j,k})) / y_width) + ((id(pressure, {i,j,max_k}) + id(pressure, {i,j,min_k})) / z_width)) - id(residual, {i,j,k})) / divide_val
          end
        end
      end end)
    res_value = :lists.sum List.flatten(dp_field)
    new_pre_field = Enum.map(Enum.to_list(0..(z_size-1)), fn(k) ->
      for j <- 0..(y_size-1) do
        for i <- 0..(x_size-1) do
          id(pressure, {i,j,k}) - id(dp_field, {i,j,k})
        end end end)
    if res_value < error_p do
      {:ok, new_pre_field}
    else
      if (ite_times+1) > max_ite_times do
        {:bad, new_pre_field}
      else
          calcModRecurse(ite_times+1, new_pre_field, residual,
            dx, dy, dx,
            x_size, y_size, z_size,
            max_ite_times, error_p, divide_val)
      end
    end
  end


  defp calcRSide {i,j,k}, {x_velocity, y_velocity, z_velocity},
    %{:dx => dx, :dy => dy, :dz => dz, :dt => dt,
    :x_size => x_size, :y_size => y_size, :z_size => z_size} do
    if 0<i && 0<j && 0<k && i<(x_size-1) && j<(y_size-1) && k<(z_size-1) do
      dudx = (id(x_velocity, {i+1,j,k}) - id(x_velocity, {i-1,j,k})) / (2 * dx)
      dvdy = (id(y_velocity, {i,j+1,k}) - id(y_velocity, {i,j-1,k})) / (2 * dy)
      dwdz = (id(z_velocity, {i,j,k+1}) - id(z_velocity, {i,j,k-1})) / (2 * dz)
      dudy = (id(x_velocity, {i,j+1,k}) - id(x_velocity, {i,j-1,k})) / (2 * dy)
      dudz = (id(x_velocity, {i,j,k+1}) - id(x_velocity, {i,j,k-1})) / (2 * dz)
      dvdx = (id(y_velocity, {i+1,j,k}) - id(y_velocity, {i-1,j,k})) / (2 * dx)
      dvdz = (id(y_velocity, {i,j,k+1}) - id(y_velocity, {i,j,k-1})) / (2 * dz)
      dwdx = (id(z_velocity, {i+1,j,k}) - id(z_velocity, {i-1,j,k})) / (2 * dx)
      dwdy = (id(z_velocity, {i,j+1,k}) - id(z_velocity, {i,j-1,k})) / (2 * dy)
      ((dudx + dvdy + dwdz) / dt) - ((dudx * dudx) + (dvdy * dvdy) + (dwdz * dwdz)) - 2*((dwdy * dvdz) + (dudz * dwdx) + (dvdx * dudy))
    else
      min_i = max 0, i-1
      max_i = min (x_size-1), i+1
      min_j = max 0, j-1
      max_j = min (y_size-1), j+1
      min_k = max 0, k-1
      max_k = min (z_size-1), k+1
      x_width = dx * (max_i - min_i)
      y_width = dy * (max_j - min_j)
      z_width = dz * (max_k - min_k)
      dudx = (id(x_velocity, {max_i,j,k}) - id(x_velocity, {min_i,j,k})) / x_width
      dvdy = (id(y_velocity, {i,max_j,k}) - id(y_velocity, {i,min_j,k})) / y_width
      dwdz = (id(z_velocity, {i,j,max_k}) - id(z_velocity, {i,j,min_k})) / z_width
      dudy = (id(x_velocity, {i,max_j,k}) - id(x_velocity, {i,min_j,k})) / y_width
      dudz = (id(x_velocity, {i,j,max_k}) - id(x_velocity, {i,j,min_k})) / z_width
      dvdx = (id(y_velocity, {max_i,j,k}) - id(y_velocity, {min_i,j,k})) / x_width
      dvdz = (id(y_velocity, {i,j,max_k}) - id(y_velocity, {i,j,min_k})) / z_width
      dwdx = (id(z_velocity, {max_i,j,k}) - id(z_velocity, {min_i,j,k})) / x_width
      dwdy = (id(z_velocity, {i,max_j,k}) - id(z_velocity, {i,min_j,k})) / y_width
      ((dudx + dvdy + dwdz) / dt) - ((dudx * dudx) + (dvdy * dvdy) + (dwdz * dwdz)) - 2*((dwdy * dvdz) + (dudz * dwdx) + (dvdx * dudy))
    end
  end


  def id enumerable, {i, j, k} do
    Enum.at(Enum.at(Enum.at(enumerable, k), j), i)
  end


  defp restrictField field, %{:x_size => x_size, :y_size => y_size, :z_size => z_size} do
    Enum.map(Enum.to_list(0..(round((z_size-1)/2)-1)), fn(k) ->
      for j <- 0..(round((y_size-1)/2)-1) do
        for i <- 0..(round((x_size-1)/2)-1) do
          (id(field, {2*i+1,2*j+1,2*k+1}) / 8 + (id(field, {2*i,2*j+1,2*k+1}) + id(field, {2*i+2,2*j+1,2*k+1}) + id(field, {2*i+1,2*j,2*k+1}) + id(field, {2*i+1,2*j+2,2*k+1}) + id(field, {2*i+1,2*j+1,2*k}) + id(field, {2*i+1,2*j+1,2*k+2})) / 16 + (id(field, {2*i,2*j,2*k+1}) + id(field, {2*i+1,2*j,2*k}) + id(field, {2*i,2*j+1,2*k}) + id(field, {2*i+2,2*j+2,2*k+1}) + id(field, {2*i+1,2*j+2,2*k+2}) + id(field, {2*i+2,2*j+1,2*k+2}) + id(field, {2*i,2*j+2,2*k+1}) + id(field, {2*i+1,2*j,2*k+2}) +  id(field, {2*i+2,2*j+1,2*k}) + id(field, {2*i+2,2*j,2*k+1}) + id(field, {2*i+1,2*j+2,2*k}) + id(field, {2*i,2*j+1,2*k+2})) / 32 + (id(field, {2*i,2*j,2*k}) + id(field, {2*i,2*j,2*k+2}) + id(field, {2*i,2*j+2,2*k}) + id(field, {2*i,2*j+2,2*k+2}) + id(field, {2*i+2,2*j,2*k}) + id(field, {2*i+2,2*j,2*k+2}) + id(field, {2*i+2,2*j+2,2*k}) + id(field, {2*i+2,2*j+2,2*k+2})) / 64)
        end
      end end)    
  end

  defp extendField field, %{:x_size => x_size, :y_size => y_size, :z_size => z_size} do
    #NOTE: use incorrect alternative method at edges and corners
    Enum.map(([1] ++ Enum.to_list(1..(z_size-2)) ++ [z_size-2]), fn(k) ->
      for j <- [1] ++ Enum.to_list(1..(y_size-2)) ++ [y_size-2] do
        for i <- [1] ++ Enum.to_list(1..(x_size-2)) ++ [z_size-2] do
          case (rem(i,2) + rem(j,2) + rem(k,2)) do
            3 ->
              id(field, {round((i-1)/2),round((j-1)/2),round((k-1)/2)})
            2 ->
              cond do
                rem(i,2) == 0 ->
                  0.5 * (id(field, {round(i/2)-1,round((j-1)/2),round((k-1)/2)}) + id(field, {round(i/2),round((j-1)/2),round((k-1)/2)}))
                rem(j,2) == 0 ->
                  0.5 * (id(field, {round((i-1)/2),round(j/2)-1,round((k-1)/2)}) + id(field, {round((i-1)/2),round(j/2),round((k-1)/2)}))
                rem(k,2) == 0 ->
                  0.5 * (id(field, {round((i-1)/2),round((j-1)/2),round(k/2)-1}) + id(field, {round((i-1)/2),round((j-1)/2),round(k/2)}))
              end
            1 ->
              cond do
                rem(i,2) == 1 ->
                  0.25 * (id(field, {round((i-1)/2),round(j/2)-1,round(k/2)-1}) + id(field, {round((i-1)/2),round(j/2)-1,round(k/2)}) + id(field, {round((i-1)/2),round(j/2),round(k/2)-1}) + id(field, {round((i-1)/2),round(j/2),round(k/2)}))
                rem(j,2) == 1 ->
                  0.25 * (id(field, {round(i/2)-1,round((j-1)/2),round(k/2)-1}) + id(field, {round(i/2)-1,round((j-1)/2),round(k/2)}) + id(field, {round(i/2),round((j-1)/2),round(k/2)-1}) + id(field, {round(i/2),round((j-1)/2),round(k/2)}))
                rem(k,2) == 1 ->
                  0.25 * (id(field, {round(i/2)-1,round(j/2)-1,round((k-1)/2)}) + id(field, {round(i/2)-1,round(j/2),round((k-1)/2)}) + id(field, {round(i/2),round(j/2)-1,round((k-1)/2)}) + id(field, {round(i/2),round(j/2),round((k-1)/2)}))
              end
            0 ->
              0.125 * (id(field, {round(i/2),round(j/2),round(k/2)}) + id(field, {round(i/2),round(j/2),round(k/2)-1}) + id(field, {round(i/2),round(j/2)-1,round(k/2)}) + id(field, {round(i/2),round(j/2)-1,round(k/2)-1}) + id(field, {round(i/2)-1,round(j/2),round(k/2)}) + id(field, {round(i/2)-1,round(j/2),round(k/2)-1}) + id(field, {round(i/2)-1,round(j/2)-1,round(k/2)}) + id(field, {round(i/2)-1,round(j/2)-1,round(k/2)-1}))
          end
        end
      end end)
  end


end # MAC.Func
