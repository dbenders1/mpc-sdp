function print_diagnostics(diagnostics)
  fprintf(" => ");
  if diagnostics.problem == 0
    cprintf([0,1,0], "feasible\n");
  else
    cprintf([1,0,0], "infeasible (problem code %i)\n",diagnostics.problem);
  end
end
