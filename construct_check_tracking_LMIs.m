function result = construct_check_tracking_LMIs(do_print,do_check,model,Q_eps,R_eps,Agrid,Bgrid,n_uxABgrid,X,Y,tol)
  result = zeros(0,2);
  n_bytes = 0;
  n_lmi = n_uxABgrid;
  n = model.n_x;
  m = model.n_u;

  for i=1:n_lmi
    A_grid = Agrid(:,:,i);
    B_grid = Bgrid(:,:,i);
    ineq = [A_grid*X+B_grid*Y+(A_grid*X+B_grid*Y)',(Q_eps*X)',(R_eps*Y)';
            Q_eps*X,-eye(n),zeros(n,m);
            R_eps*Y,zeros(m,n),-eye(m)];
    if ~do_check
      result = [result;ineq<=0];
      if do_print && rem(i,200) == 0
        fprintf(repmat('\b',1,n_bytes));
        n_bytes = fprintf("Constructed tracking LMI %i/%i\n",i,n_lmi);
      end
    else
      max_eig(i) = max(eig(ineq));
      if (max_eig(i) > tol)
        result = [result;[i,max_eig(i)]];
        if do_print
          fprintf("Maximum eigenvalue of tracking LMI %i > tol!\n",i);
        end
      end
      if do_print && rem(i,200) == 0
        fprintf(repmat('\b',1,n_bytes));
        n_bytes = fprintf("Checked tracking LMI %i/%i\n",i,n_lmi);
      end
    end
  end

  if ~do_check
    if do_print
      fprintf(repmat('\b',1,n_bytes));
      fprintf("Constructed tracking LMI %i/%i\n",i,n_lmi);
    end
  else
    if do_print
      fprintf(repmat('\b',1,n_bytes));
      fprintf("Checked tracking LMI %i/%i\n",i,n_lmi);
    end
  end
end
