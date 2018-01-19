defmodule CalcPServerTest do
  use ExUnit.Case
  doctest CalcPServer
  @tag timeout: 3600000

  test "calculate next step's pressure" do
    v_field = List.duplicate([0] ++ List.duplicate(1, 399) ++ [2], 201)
    p_field = List.duplicate([0] ++ List.duplicate(1, 399) ++ [2], 201)
    bc_field = for j <- 0..200 do
      for i <- 0..400 do
        if 49 == i || 49 == j do
          10
        else
          false
        end
      end end
    CalcPServer.genCalcServer %{:max_ite_times => 100,
                                :error_p => 0.0001,
                                :omega => 1,
                                :max_res_ratio => 0.5}
    {status, result} = CalcPServer.calcPre {v_field, v_field}, p_field, bc_field,
      %{:x_size => 401,
        :y_size => 201,
        :dx => 0.1,
        :dy => 0.1,
        :dt => 0.01,
        :Re => 70}
    # IO.inspect result
    IO.inspect status
    
  end
end
