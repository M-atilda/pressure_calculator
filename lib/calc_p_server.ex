defmodule CalcPServer do
  @moduledoc """
  Documentation for CalcPServer.
  """
  @doc """
  Hello world.

  ## Examples

      iex> CalcPServer.hello
      :world

  """
  @name p_server_name :g_p_calc_server
  import MAC.Func

  def calcPre velocitys_field, pressure_field, bc_field, information do
    server_name = :global.whereis_name(@name)
    send server_name, {:calc, velocitys_field, pressure_field, bc_field, information, self}
    receive do
      {simbol, result, ^server_name} ->
        case simbol do
          :error ->
            raise RuntimeError
          _ ->
            #TODO: when receive :bad(calculation hasnot finished in max iteration times)
            result
        end
    end
  end

  def genCalcServer calc_info do
    pid = spawn(__MODULE__, :calc_server, [calc_info])
    :global.register_name(@name, pid)
  end

  def calc_server calc_info do
    receive do
      {:calc, velocitys_field, pressure_field, bc_field, information, client} ->
        {status, result} = derivePre velocitys_field, pressure, bc_field, information, calc_info
        send client, {status, result, self}
    end
    calc_server calc_info
  end

end # CalcPServer
