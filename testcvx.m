cvx_begin sdp quiet
     %   cvx_solver Mosek
                  
            variable powerAllocateXstar(2 , 2,3) hermitian semidefinite
        
      %  variable totalPower nonnegative
        variable gamma2p nonnegative
        expression littleA
        minimize gamma2p
        subject to
            for iUser = 1:3
                real(  trace(eye(2) * powerAllocateXstar(:,:,iUser) ) )- 1 - gamma2p <= 0;
                littleA  = littleA + real(  trace(eye(2) * powerAllocateXstar(:,:,iUser) ));
            end
            littleA = 1;
   cvx_end