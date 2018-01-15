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
  def hello do
    :world
  end
  @name :g_p_calc_server
  import MAC.Func
  

  #TODO: make server for each calculation (and give them original name)
  def calcPre velocitys_field, pressure_field, bc_field, information do
    server = :global.whereis_name(@name)
    send server, {:calc, velocitys_field, pressure_field, bc_field, information, self}
    receive do
      {simbol, result, ^server} ->
        case simbol do
          :error ->
            raise RuntimeError
          status ->
            #TODO: when receive :bad(calculation hasnot finished in max iteration times)
            {status, result}
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
        {status, result} = derivePre velocitys_field, pressure_field, bc_field, information, calc_info
        send client, {status, result, self}
    end
    calc_server calc_info
  end

end # CalcPServer
