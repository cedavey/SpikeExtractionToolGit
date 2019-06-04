function bool = isint(x)
  test = @(x) round(x) == x;
  bool = all(arrayfun(test,x));
end