defmodule CalcPServerTest do
  use ExUnit.Case
  doctest CalcPServer
  @tag timeout: 3600000

  test "calculate next step's pressure" do
    line = List.to_tuple [0.0|List.duplicate(1.0, 399) ++ [2]]
    v_field = Tuple.duplicate(line, 201)
    p_field = Tuple.duplicate(line, 201)
    bc_field = for j <- 0..200 do
      for i <- 0..400 do
        if 49 == i || 49 == j do
          10
        else
          false
        end
      end
      |> List.to_tuple
    end
    |> List.to_tuple
    CalcPServer.genCalcServer "test", %{:max_ite_times => 1000,
                                        :error_p => 0.0001,
                                        :omega => 1.0,
                                        :max_res_ratio => 1}
    {status, result} = CalcPServer.calcPre {v_field, v_field}, p_field, bc_field,
      %{:x_size => 401,
        :y_size => 201,
        :dx => 0.1,
        :dy => 0.1,
        :dt => 0.01,
        :Re => 70}, "test"
    IO.inspect result
    IO.inspect status
    
  end
end
