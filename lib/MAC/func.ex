#file   func.ex
#author mi-na
#date   18/01/14
#brief  calculate pressure field
#       right side is derived using central difference
#       in recursion, it applies 2 level multi-grid scheme(v cycle) and SOR method
defmodule MAC.Func do


  def derivePre velocitys_field, pressure, bc_field,
    %{:dx => dx,
      :dy => dy,
      :dz => dz,
      :x_size => x_size,
      :y_size => y_size,
      :z_size => z_size}=information,
    %{:max_ite_times => max_ite_times,
      :error_p => error_p,
      :omega => omega} do
    right_side = for k <- 0..(z_size-1) do
      for j <- 0..(y_size-1) do
        for i <- 0..(x_size-1) do
          if 0<i && 0<j && 0<k && i<(x_size-1) && j<(y_size-1) && k<(z_size-1) do
            calcRSide {i,j,k}, velocitys_field, information
          else
            0
          end
        end
      end
    end
    IO.inspect right_side
    {status, new_pressure, residual_field} = derivePreRecurse(0, right_side, pressure, nil, bc_field,
      information,
      max_ite_times, error_p, omega)
    {status, new_pressure}
  end
  def derivePreRecurse ite_times, right_side, pressure, residual, bc_field,
    %{:dx => dx,
      :dy => dy,
      :dz => dz,
      :x_size => x_size,
      :y_size => y_size,
      :z_size => z_size}=information,
    max_ite_times, error_p, omega do
    #NOTE: these calculation information calues are arranged in derivePre for multi-grid scheme
    if ite_times > max_ite_times do
      {:bad, pressure, residual}
    else
      {has_calced, {new_pressure, new_residual, residual_value}} = try do
                                                                     {true, derivePreStep(pressure, right_side, bc_field, dx,dy,dx,x_size,y_size,z_size, omega)}
                                                                   rescue
                                                                     err ->
                                                                       IO.puts "[Error] #{inspect err}"
                                                                     {false, {pressure, residual, 0}}
                                                                   end
      if has_calced do
        if residual_value < error_p do
          {:ok, new_pressure, new_residual}
        else
          IO.inspect new_pressure
          derivePreRecurse(ite_times+1, right_side, new_pressure, new_residual, bc_field,
            information,
            max_ite_times, error_p, omega)
        end
      else
        {:error, pressure, residual}
      end
    end
  end

  def derivePreStep pressure, right_side, bc_field,
    dx, dy, dz,
    x_size, y_size, z_size,
    omega do
    divide_val = 2*(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz))
    new_pressure = for k <- 0..(z_size-1) do
      for j <- 0..(y_size-1) do
        for i <- 0..(x_size-1) do
          if id(bc_field, {i,j,k}) == nil do
            if 0<i && 0<j && 0<k && i<(x_size-1) && j<(y_size-1) && k<(z_size-1) do
              dp = ((((id(pressure, {i+1,j,k}) + id(pressure, {i-1,j,k})) / (2*dx))
                + ((id(pressure, {i,j+1,k}) + id(pressure, {i,j-1,k})) / (2*dy))
                + ((id(pressure, {i,j,k+1}) + id(pressure, {i,j,k-1})) / (2*dz))
                - id(right_side, {i,j,k})) /
                divide_val) -
              id(pressure, {i,j,k})
              {id(pressure, {i,j,k}) + omega*dp, dp*dp}
            else
              {id(pressure, {i,j,k}), 0}
            end
          else
            {id(bc_field, {i,j,k}), 0}
          end
        end
      end
    end
    residual = Enum.reduce List.flatten(new_pressure), 0, fn({_p, dr}, acm) -> acm+dr end
    {for k <- 0..(z_size-1) do
         for j <- 0..(y_size-1) do
           for i <- 0..(x_size-1) do
             {p, _dr} = id(new_pressure, {i,j,k})
             p
           end
         end
     end,
     for k <- 0..(z_size-1) do
         for j <- 0..(y_size-1) do
           for i <- 0..(x_size-1) do
             {_p, dr} = id(new_pressure, {i,j,k})
             dr
           end
         end
     end,
     residual}
  end


  def calcRSide {i,j,k}, {x_velocity, y_velocity, z_velocity}, %{:dx => dx, :dy => dy, :dz => dz, :dt => dt} do
    dudx = (id(x_velocity, {i+1,j,k}) - id(x_velocity, {i-1,j,k})) / (2 * dx)
    dvdy = (id(y_velocity, {i,j+1,k}) - id(y_velocity, {i,j-1,k})) / (2 * dx)
    dwdz = (id(z_velocity, {i,j,k+1}) - id(z_velocity, {i,j,k-1})) / (2 * dx)
    dudy = (id(x_velocity, {i,j+1,k}) - id(x_velocity, {i,j-1,k})) / (2 * dx)
    dudz = (id(x_velocity, {i,j,k+1}) - id(x_velocity, {i,j,k-1})) / (2 * dx)
    dvdx = (id(y_velocity, {i+1,j,k}) - id(y_velocity, {i-1,j,k})) / (2 * dx)
    dvdz = (id(y_velocity, {i,j,k+1}) - id(y_velocity, {i,j,k-1})) / (2 * dx)
    dwdx = (id(z_velocity, {i+1,j,k}) - id(z_velocity, {i-1,j,k})) / (2 * dx)
    dwdy = (id(z_velocity, {i,j+1,k}) - id(z_velocity, {i,j-1,k})) / (2 * dx)
    ((dudx + dvdy + dwdz) / dt) -
    ((dudx * dudx) + (dvdy * dvdy) + (dwdz * dwdz)) -
    2*((dwdy * dvdz) + (dudz * dwdx) + (dvdx * dudy))
  end

  def id enumerable, {i, j, k} do
    Enum.at(Enum.at(Enum.at(enumerable, k), j), i)
  end


end # MAC.Func
