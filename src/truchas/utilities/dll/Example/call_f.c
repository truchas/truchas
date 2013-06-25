float call_f_(float (**addr)(), float *x)
{
  return (*addr)(x);
}
