clear;clc
cvx_begin sdp quiet
        cvx_solver Mosek
                  
            variable powerAllocateXstar(8 , 8,  3) hermitian semidefinite
       
      %  variable totalPower nonnegative
        variable gamma2p nonnegative
        %minimize gamma2p
        subject to
            
                for iUser = 1:3
            ( trace(powerAllocateXstar(:,:,iUser)))   == 1 ;
                end
   cvx_end