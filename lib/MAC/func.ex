#file   func.ex
#author mi-na
#date   18/01/14
#brief  calculate pressure field
#       right side is derived using central difference
#       in recursion, it applies 2 level multi-grid scheme(v cycle) and SOR method
defmodule MAC.Func do


  def derivePre velocitys_field, pressure, bc_field,
    %{:x_size => x_size,
      :y_size => y_size}=information,
    %{:max_ite_times => max_ite_times,
      :error_p => error_p,
      :omega => omega,
      :max_res_ratio => max_res_ratio} do
    x_half_size = round((x_size-1) / 2)
    y_half_size = round((y_size-1) / 2)
    left_up_dp_task = Task.async(fn -> calcRSidePartially {0..x_half_size, 0..y_half_size}, velocitys_field, information end)
    left_down_dp_task = Task.async(fn -> calcRSidePartially {0..x_half_size, (y_half_size+1)..(y_size-1)}, velocitys_field, information end)
    right_up_dp_task = Task.async(fn -> calcRSidePartially {(x_half_size+1)..(x_size-1), 0..y_half_size}, velocitys_field, information end)
    right_down_dp_task = Task.async(fn -> calcRSidePartially {(x_half_size+1)..(x_size-1), (y_half_size+1)..(y_size-1)}, velocitys_field, information end)
    left_up_dp = Task.await(left_up_dp_task)
    left_down_dp = Task.await(left_down_dp_task)
    right_up_dp = Task.await(right_up_dp_task)
    right_down_dp = Task.await(right_down_dp_task)
    right_side = Enum.map(Enum.to_list(0..(y_size-1)), fn(j) ->
      for i <- 0..(x_size-1) do
        cond do
          i<(x_half_size+1) && j<(y_half_size+1) ->
            id(left_up_dp, {i,j})
          i<(x_half_size+1) && j>y_half_size ->
            id(left_down_dp, {i,j-y_half_size-1})
          i>x_half_size && j<(y_half_size+1) ->
            id(right_up_dp, {i-x_half_size-1,j})
          i>x_half_size && j>y_half_size ->
            id(right_down_dp, {i-x_half_size-1,j-y_half_size-1})
        end
      end
    end)
    derivePreRecurse(0, right_side, pressure, bc_field,
      information,
      max_ite_times, error_p, omega,
      max_res_ratio, 1, 1)
  end
  def derivePreRecurse ite_times, right_side, pressure, bc_field,
    %{:dx => dx,
      :dy => dy,
      :x_size => x_size,
      :y_size => y_size}=information,
    max_ite_times, error_p, omega,
    max_res_ratio, res_value, res_ratio do
    if ite_times > max_ite_times do
      {:bad, pressure}
    else
      {has_calced, {new_pressure, new_res_field, new_res_value}} = try do
                                                                     {true, derivePreStep(pressure, right_side, bc_field, dx,dy,x_size,y_size, omega)}
                                                                   rescue
                                                                     err ->
                                                                       IO.puts "[Error] #{inspect err}"
                                                                     {false, {pressure, nil, 0}}
                                                                   end
      if has_calced do
        try do
          IO.puts "(P) #{round(:math.log(new_res_value) / :math.log(error_p) * 100)}% #{inspect DateTime.utc_now}"
        rescue
          _ -> nil
        end
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
            modified_new_pressure = Enum.map(Enum.to_list(0..(y_size-1)), fn(j) ->
              for i <- 0..(x_size-1) do
                id(new_pressure, {i,j}) + id(dp_field, {i,j})
              end end)
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
    dx, dy,
    x_size, y_size,
    omega do
    divide_val = 2*(1/(dx*dx) + 1/(dy*dy))


    # divide to 4 spaces for pararell processing
    x_half_size = round((x_size-1) / 2)
    y_half_size = round((y_size-1) / 2)
    left_up_dp_task = Task.async(fn -> deriveDPPartially pressure, right_side, bc_field, divide_val, {dx,dy}, {x_size,y_size}, {0..x_half_size, 0..y_half_size} end)
    left_down_dp_task = Task.async(fn -> deriveDPPartially pressure, right_side, bc_field, divide_val, {dx,dy}, {x_size,y_size}, {0..x_half_size, (y_half_size+1)..(y_size-1)} end)
    right_up_dp_task = Task.async(fn -> deriveDPPartially pressure, right_side, bc_field, divide_val, {dx,dy}, {x_size,y_size}, {(x_half_size+1)..(x_size-1), 0..y_half_size} end)
    right_down_dp_task = Task.async(fn -> deriveDPPartially pressure, right_side, bc_field, divide_val, {dx,dy}, {x_size,y_size}, {(x_half_size+1)..(x_size-1), (y_half_size+1)..(y_size-1)} end)
    left_up_dp = Task.await(left_up_dp_task)
    left_down_dp = Task.await(left_down_dp_task)
    right_up_dp = Task.await(right_up_dp_task)
    right_down_dp = Task.await(right_down_dp_task)
    dp_field = Enum.map(Enum.to_list(0..(y_size-1)), fn(j) ->
      for i <- 0..(x_size-1) do
        cond do
          i<(x_half_size+1) && j<(y_half_size+1) ->
            id(left_up_dp, {i,j})
          i<(x_half_size+1) && j>y_half_size ->
            id(left_down_dp, {i,j-y_half_size-1})
          i>x_half_size && j<(y_half_size+1) ->
            id(right_up_dp, {i-x_half_size-1,j})
          i>x_half_size && j>y_half_size ->
            id(right_down_dp, {i-x_half_size-1,j-y_half_size-1})
        end
      end
    end)


  {Enum.map(Enum.to_list(0..(y_size-1)), fn(j) ->
      for i <- 0..(x_size-1) do
        if id(bc_field, {i,j}) == nil do
          id(pressure, {i,j}) + omega * id(dp_field, {i,j})
        else
          id(bc_field, {i,j})
        end
      end
    end),
   dp_field,
   List.flatten(dp_field)
   |> Enum.map(&(&1*&1))
   |> :lists.sum}
  end
  defp deriveDPPartially pressure, right_side, bc_field, divide_val, {dx,dy}, {x_size,y_size}, {x_range,y_range} do
    Enum.map(y_range, fn(j) ->
      for i <- x_range do
        if id(bc_field, {i,j}) == nil do
          if 0<i && 0<j && i<(x_size-1) && j<(y_size-1) do
            (((id(pressure, {i+1,j}) + id(pressure, {i-1,j})) / (2*dx)) + ((id(pressure, {i,j+1}) + id(pressure, {i,j-1})) / (2*dy)) - id(right_side, {i,j})) / divide_val - id(pressure, {i,j})
          else
              min_i = max 0, i-1
              max_i = min (x_size-1), i+1
              min_j = max 0, j-1
              max_j = min (y_size-1), j+1
              x_width = dx * (max_i - min_i)
              y_width = dy * (max_j - min_j)
              (((id(pressure, {max_i,j}) + id(pressure, {min_i,j})) / x_width) + ((id(pressure, {i,max_j}) + id(pressure, {i,min_j})) / y_width) - id(right_side, {i,j})) / divide_val - id(pressure, {i,j})
            end
          else
            0
          end
        end
    end)
  end
  
  def calcModField residual,
    %{:dx => small_dx,
      :dy => small_dy,
      :x_size => small_x_size,
      :y_size => small_y_size}, max_ite_times, error_p do
    IO.puts "[Info] start calcModField #{inspect DateTime.utc_now}"    
    dx = 2 * small_dx
    dy = 2 * small_dy
    x_size = round((small_x_size - 1) / 2)
    y_size = round((small_y_size - 1) / 2)
    divide_val = 2*(1/(dx*dx) + 1/(dy*dy))
    zero_pressure = Enum.map(Enum.to_list(0..(y_size-1)), fn(_j) ->
      for _i <- 0..(x_size-1) do
        0
      end end)
    {status, new_pre_field} = try do
                               calcModRecurse(0, zero_pressure, residual,
                                 dx, dy,
                                 x_size, y_size,
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
    dx, dy,
    x_size, y_size,
    max_ite_times, error_p, divide_val do
    dp_field = Enum.map(Enum.to_list(0..(y_size-1)), fn(j) ->
      for i <- 0..(x_size-1) do
        if 0<i && 0<j && i<(x_size-1) && j<(y_size-1) do
          ((((id(pressure, {i+1,j}) + id(pressure, {i-1,j})) / (2*dx)) + ((id(pressure, {i,j+1}) + id(pressure, {i,j-1})) / (2*dy))) - id(residual, {i,j})) / divide_val
        else
          min_i = max 0, i-1
          max_i = min (x_size-1), i+1
          min_j = max 0, j-1
          max_j = min (y_size-1), j+1
          x_width = dx * (max_i - min_i)
          y_width = dy * (max_j - min_j)
          ((((id(pressure, {max_i,j}) + id(pressure, {min_i,j})) / x_width) + ((id(pressure, {i,max_j}) + id(pressure, {i,min_j})) / y_width) - id(residual, {i,j}))) / divide_val
        end
      end
    end)
    res_value = :lists.sum List.flatten(dp_field)
    new_pre_field = Enum.map(Enum.to_list(0..(y_size-1)), fn(j) ->
      for i <- 0..(x_size-1) do
        id(pressure, {i,j}) - id(dp_field, {i,j})
      end end)
    if res_value < error_p do
      {:ok, new_pre_field}
    else
      if (ite_times+1) > max_ite_times do
        {:bad, new_pre_field}
      else
          calcModRecurse(ite_times+1, new_pre_field, residual,
            dx, dy,
            x_size, y_size,
            max_ite_times, error_p, divide_val)
      end
    end
  end
  

  defp calcRSide {i,j}, {x_velocity, y_velocity},
    %{:dx => dx, :dy => dy, :dt => dt,
      :x_size => x_size, :y_size => y_size} do
    if 0<i && 0<j && i<(x_size-1) && j<(y_size-1) do
      dudx = (id(x_velocity, {i+1,j}) - id(x_velocity, {i-1,j})) / (2 * dx)
      dvdy = (id(y_velocity, {i,j+1}) - id(y_velocity, {i,j-1})) / (2 * dy)
      dudy = (id(x_velocity, {i,j+1}) - id(x_velocity, {i,j-1})) / (2 * dy)
      dvdx = (id(y_velocity, {i+1,j}) - id(y_velocity, {i-1,j})) / (2 * dx)
      ((dudx + dvdy) / dt) - ((dudx * dudx) + (dvdy * dvdy)) - 2*(dvdx * dudy)
    else
      min_i = max 0, i-1
      max_i = min (x_size-1), i+1
      min_j = max 0, j-1
      max_j = min (y_size-1), j+1
      x_width = dx * (max_i - min_i)
      y_width = dy * (max_j - min_j)
      dudx = (id(x_velocity, {max_i,j}) - id(x_velocity, {min_i,j})) / x_width
      dvdy = (id(y_velocity, {i,max_j}) - id(y_velocity, {i,min_j})) / y_width
      dudy = (id(x_velocity, {i,max_j}) - id(x_velocity, {i,min_j})) / y_width
      dvdx = (id(y_velocity, {max_i,j}) - id(y_velocity, {min_i,j})) / x_width
      ((dudx + dvdy) / dt) - ((dudx * dudx) + (dvdy * dvdy)) - 2*(dvdx * dudy)
    end
  end
  defp calcRSidePartially {x_range, y_range}, velocitys_field, information do
    for j <- y_range do
      for i <- x_range do
        calcRSide {i,j}, velocitys_field, information
      end
    end
  end

  def id enumerable, {i, j} do
    Enum.at(Enum.at(enumerable, j), i)
  end
  

  defp restrictField field, %{:x_size => x_size, :y_size => y_size} do
    Enum.map(Enum.to_list(0..(round((y_size-1)/2)-1)), fn(j) ->
      for i <- 0..(round((x_size-1)/2)-1) do
        id(field, {2*i+1,2*j+1}) / 4 + (id(field, {2*i,2*j+1}) + id(field, {2*i+2,2*j+1}) + id(field, {2*i+1,2*j}) + id(field, {2*i+1,2*j+2})) / 8 + (id(field, {2*i,2*j}) + id(field, {2*i,2*j+2}) + id(field, {2*i+2,2*j}) + id(field, {2*i+2,2*j+2})) / 16
      end
    end)    
  end

  defp extendField field, %{:x_size => x_size, :y_size => y_size} do
    #NOTE: use incorrect alternative method at edges and corners
    Enum.map(([1] ++ Enum.to_list(1..(y_size-2)) ++ [y_size-2]), fn(j) ->
      for i <- [1] ++ Enum.to_list(1..(x_size-2)) ++ [x_size-2] do
        case (rem(i,2) + rem(j,2)) do
          2 ->
            id(field, {round((i-1)/2),round((j-1)/2),})
          1 ->
            cond do
              rem(i,2) == 0 ->
                0.5 * (id(field, {round(i/2)-1,round((j-1)/2)}) + id(field, {round(i/2),round((j-1)/2)}))
              rem(j,2) == 0 ->
                0.5 * (id(field, {round((i-1)/2),round(j/2)-1}) + id(field, {round((i-1)/2),round(j/2)}))
            end
          0 ->
            0.25 * (id(field, {round(i/2),round(j/2)}) + id(field, {round(i/2),round(j/2)-1}) + id(field, {round(i/2)-1,round(j/2)}) + id(field, {round(i/2)-1,round(j/2)-1}))
        end
      end
    end)
  end


end # MAC.Func
