function result = construct_check_sys_constraint_LMIs(do_print,do_check,n_s,L_s,eps_s,X,Y,tol)
  result = zeros(0,2);
  n_bytes = 0;
  n_lmi = n_s;

  % Initialize result with no LMIs being violated, fill below if violated
  if do_check
    result = zeros(0,2);
  end

  % Check computation of eps_s, related to minimizing size of RPI set
  i = 0;
  for k = 1:n_s
    i = i + 1;
    ineq = [eps_s(k),L_s(k,:)*[Y;X];
            (L_s(k,:)*[Y;X])',X];

    if ~do_check
      result = [result;ineq>=0];
      if do_print && rem(i,200) == 0
        fprintf(repmat('\b',1,n_bytes));
        n_bytes = fprintf("Constructed system constraint LMI %i/%i\n",i,n_lmi);
      end
    else
      min_eig(i) = min(eig(ineq));
      if (min_eig(i) < -tol)
        result = [result;[i,min_eig(i)]];
        if do_print
          fprintf("Minimum eigenvalue of system constraint LMI %i < -tol!\n",i);
        end
      end
      if do_print && rem(i,200) == 0
        fprintf(repmat('\b',1,n_bytes));
        n_bytes = fprintf("Checked system constraint LMI %i/%i\n",i,n_lmi);
      end
    end
  end

  if ~do_check
    if do_print
      fprintf(repmat('\b',1,n_bytes));
      fprintf("Constructed system constraint LMI %i/%i\n",i,n_lmi);
    end
  else
    if do_print
      fprintf(repmat('\b',1,n_bytes));
      fprintf("Checked system constraint LMI %i/%i\n",i,n_lmi);
    end
  end
end
