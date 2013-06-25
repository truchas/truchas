#     Author: John A. Turner
#             e-mail: john.turner@pobox.com
#     $Id: summary.awk,v 1.2 2002/12/04 00:57:57 turner Exp $
#
#     awk script to extract summary from UbikTest output.

{
  if ($1 == "Max." && $5 == "solution")
	error = $7
  else if ($1 == "Iterations")
    iters = $3
  else if ($1 == "CPU")
    cpu = $7
  else if ($1 == "Storage")
  {
    if ($4 == "ELLPACK-ITPACK")
      storage = "ELL "
    else if ($4 == "Coordinate")
      storage = "COO "
    else if ($4 == "Row")
      storage = "RSS "
    else if ($4 == "Column")
      storage = "CSS "
    else
      storage = $4
  }
  else if ($1 == "Method")
  {
    if ($4 == "GMRES")
      method = $3 " " $4 " " $5
    else
      method = $3
  }
  else if ($1 == "Preconditioner:")
  {
    if ($2 == "m-step")
      precond = $3
    else
      precond = $2
  }
  else if ($1 == "System")
    scaling = " / " $3
  else if ($1 == "Jacobi")
    scaling = scaling ", diagonal scaled"
  else if ($1 == "Relaxation")
    omega = $3
  else if ($1 == "Number")
  {
	if ($3 == "unknowns")
      unknowns = $5
	else if ($3 == "steps:")
	  steps = $4
  }
  else if ($1 == "Block")
    blocksize = $4
  else if ($1 == "Krylov")
    subspace = $4
  else if ($1 == "Convergence" && $3 == "not")
    error_msg = $0
}
END {

  if (method == "")
    exit
  else if (method == "LU")
    method = method scaling
  else if (method == "Jacobi" || method == "SOR" || method == "SSOR")
    method = method "(" omega ")" scaling
  else
  {
    if (method == "GMRES")
	  method = method "(" subspace ")"

    if (precond == "Jacobi" || precond == "SSOR") 
      precond = precond "(" steps "/" omega ")"
    else if (precond == "IC" || precond == "ILU")
    {
      precond = precond "(" omega ")"
      method = method " / " precond scaling
    }
  }

  if (cpuout == "yes")
    printf "%8.3f %s %s\n", cpu, " |  " storage " /", method
  else if (perfout == "yes")
    printf "    %5s (%4s)       %5s     %8.3f   %s\n", unknowns, blocksize, iters, cpu, error
  else
  {
	printf "%5s  %s %s %s\n", iters, error, " |  " storage " /", method
    if (error_msg != "") printf "   (%s)\n", error_msg
  }
}
