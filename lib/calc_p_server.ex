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
  def hello, do: :world
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
    IO.puts "[Info] start calc_P_server <#{inspect pid}>"
  end
  def genCalcServer do
    genCalcServer %{:max_ite_times => 100,
                    :error_p => 0.0001,
                    :omega => 1,
                    :max_res_ratio => 0.5}
  end

  def calc_server calc_info do
    receive do
      {:calc, velocitys_field, pressure_field, bc_field, information, client} ->
        IO.puts "[Info] pressure calculation <#{inspect client}> #{inspect DateTime.utc_now}"
        {status, result} = derivePre velocitys_field, pressure_field, bc_field, information, calc_info
        IO.puts "[Info] calculation finished <#{inspect client}> #{inspect DateTime.utc_now}"
        send client, {status, result, self}
    end
    calc_server calc_info
  end

end # CalcPServer
