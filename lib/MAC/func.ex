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
      :error_p => error_p}=calc_info do
    right_side = for k <- 0..z_size do
      for j <- 0..y_size do
        for i <- 0..x_size do
          if 0<i && 0<j && 0<k && i<x_size && j<y_size && k<z_size do
            calcRSide {i,j,k}, velocitys_field, information
          else
            0
          end
        end
      end
    end

    derive_pressure_recursive_fn = fn(current_pressure, ite_times) ->
      if ite_times > max_ite_times do
        {:bad, current_pressure}
      else
        {has_calced, next_pressure, residual} = try do
                                                  {true, derivePreStep(current_pressure, right_side, bc_field, dx,dy,dx,x_size,y_size,z_size, calc_info)}
                                                rescue
                                                  _ -> {false, {current_pressure, 0}}
                                                end
        if has_calced do
          if residual < error_p do
            {:ok, next_pressure}
          else
            derive_pressure_recursive_fn next_pressure, ite_times+1
          end
        else
          {:error, current_pressure}
        end
      end
    end
    derive_pressure_recursive_fn pressure, 0
  end

  def derivePreStep pressure, right_side, bc_field,
    dx, dy, dz,
    x_size, y_size, z_size,
    %{:omega => omega} do
    divide_val = 2*(1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz))
    new_pressure = for k <- 0..z_size do
      for j <- 0..y_size do
        for i <- 0..x_size do
          if id(bc_field, {i,j,k}) == nil do
            if 0<i && 0<j && 0<k && i<x_size && j<y_size && k<z_size do
              dp = (((id(pressure, {i+1,j,k}) + id(pressure, {i-1,j,k})) / (2*dx))
                + ((id(pressure, {i,j+1,k}) + id(pressure, {i,j-1,k})) / (2*dy))
                + ((id(pressure, {i,j,k+1}) + id(pressure, {i,j,k-1})) / (2*dz))
                + id(right_side, {i,j,k})) /
              divide_val
              {id(pressure, {i,j,k}) + omega*dp, dp*dp}
            else
              id(pressure, {i,j,k})
            end
          else
            id(bc_field, {i,j,k})
          end
        end
      end
    end
    residual = Enum.reduce List.flatten(new_pressure), 0, fn({_p, dr}, acm) -> acm+dr end
    {for k <- 0..z_size do
         for j <- 0..y_size do
           for i <- 0..x_size do
             id(new_pressure, {i,j,k})
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
