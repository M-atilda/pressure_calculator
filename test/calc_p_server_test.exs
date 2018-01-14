defmodule CalcPServerTest do
  use ExUnit.Case
  doctest CalcPServer

  test "greets the world" do
    assert CalcPServer.hello() == :world
  end
end
