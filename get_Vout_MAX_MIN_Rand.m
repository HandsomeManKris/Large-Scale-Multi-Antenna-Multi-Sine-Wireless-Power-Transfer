function [sprecoder,VoutUser,Vout] =  get_Vout_MAX_MIN_Rand(ChannelInfo,TransceiverInfo,InitialM)






T = TransceiverInfo.T;
GainFre = ChannelInfo.subbandChannelGainFre;
Mt = TransceiverInfo.Mt;
Tolerance = TransceiverInfo.Tolerance;
MrPower = TransceiverInfo.MrPower;
K =  TransceiverInfo.K;
b2 = TransceiverInfo.b2;
b4 = TransceiverInfo.b4;
Mdiag = InitialM.Mdiag;
subbandNumber = ChannelInfo.subbandNumber;
powerAllocateXlast = InitialM.powerAllocateX;
tk = InitialM.tk;
tkstar = InitialM.tkstar;
Loop = 1;
wq = 1;

A0 = diag([0.5 ones(1,subbandNumber-1)] * (-3 * b4));

while(Loop) 
    %update cqbar Cq1,Aq1
    gammaGroup = zeros(K,T);
    cqbar = zeros(K,1);
    Cq1 = zeros(subbandNumber * Mt,subbandNumber *  Mt,K);
    Aq1 = zeros(subbandNumber * Mt,subbandNumber *  Mt,K);
    for iUser  = 1:K
        cqbar(iUser) = -tk(:,iUser)' * A0 * tk(:,iUser);
         Cq1(:,:,iUser) =  -1/2 * (b2 + 3 * b4 * tk(1,iUser)) * Mdiag(:,:,1,iUser);
         for iSubbandNumber = 2:subbandNumber
             Cq1(:,:,iUser) = Cq1(:,:,iUser) - 3* b4 * tkstar(iSubbandNumber-1,iUser) * Mdiag(:,:,iSubbandNumber,iUser);
         end
          Aq1(:,:,iUser) = Cq1(:,:,iUser) + Cq1(:,:,iUser)';
    end
       
     
    cvx_begin sdp quiet
        cvx_solver Mosek
        
        variable  powerAllocateX(subbandNumber * Mt,subbandNumber * Mt)   semidefinite
        %variable  powerAllocateX(subbandNumber * Mt,subbandNumber * Mt) hermitian 
        variable gamma0  nonnegative
        minimize -gamma0
        subject to 
            for iUser = 1:K
               real(trace(Aq1(:,:,iUser) * powerAllocateX) + cqbar(iUser)) + gamma0 <= 0;
                
            end
           % real(trace(Aq1(:,:,1) * powerAllocateX) + cqbar(1)) + gamma0 <= 0;
           % real(trace(Aq1(:,:,2) * powerAllocateX) + cqbar(2)) + gamma0 <= 0;
          %  real(trace(Aq1(:,:,3) * powerAllocateX) + cqbar(3) + gamma0 <= 0;
            trace( powerAllocateX)<= MrPower;
    cvx_end
    
   % for iUser = 1:K
   %     gammaGroup(iUser)  =   trace(Aq1(:,:,iUser) * powerAllocateX) + cqbar(iUser);
   % end
   
   % [~,userIndex] = max(gammaGroup);
   %create T vectors
   a = zeros(Mt * subbandNumber,1);
   vect = zeros(Mt * subbandNumber, T);
   vectT = zeros(Mt * subbandNumber, T);
   powerAllocateXGroup = zeros(Mt * subbandNumber,Mt * subbandNumber,T);
   negGamma = zeros(T,1);
   [VX,DX] = eig( powerAllocateX);
   for iT = 1:T
       a = rand(Mt * subbandNumber,1);
       vect(:,iT) = cos(a * 2 * pi) + sin(a * 2 * pi) * 1i;
       vectT(:,iT) = VX * DX^(1/2) * vect(:,iT);
       powerAllocateXGroup(:,:,iT) = vectT(:,iT) * vectT(:,iT)';
       for iUser =  1:K
           gammaGroup(iUser,iT)  =   trace(Aq1(:,:,iUser) * powerAllocateXGroup(:,:,iT)) + cqbar(iUser);
       end
       negGamma(iT) = max(gammaGroup(:,iT));
   end
   [~,minIndex] = min(negGamma);
   
    xstar =  VX * DX^(1/2) *   vect(:,minIndex);
     powerAllocateXstar = xstar * xstar';
    %Update tk tkstar
    tk = zeros(subbandNumber,K);
    tkstar = zeros(subbandNumber-1,K);
    for iUser = 1:K
        for iSubbandNumber = 1:subbandNumber
            %Mdiag(:,:, iSubbandNumber,iUser) = diag(diag(Mpp(:,:,iUser),iSubbandNumber-1),iSubbandNumber-1);
            tk(iSubbandNumber,iUser) = trace(Mdiag(:,:, iSubbandNumber,iUser) * powerAllocateXstar);
            if iSubbandNumber ~= 1
                tkstar(iSubbandNumber-1,iUser) = trace( Mdiag(:,:,iSubbandNumber,iUser)' * powerAllocateXstar);
            end
        end
    
    end
      norm((powerAllocateXstar-powerAllocateXlast),'fro')/norm(powerAllocateXstar,'fro')
     if (norm((powerAllocateXstar-powerAllocateXlast),'fro')/norm(powerAllocateXstar,'fro')) <= Tolerance/10 || Loop>=20
        
        break;
    else
         Loop = Loop+1;
        
        powerAllocateXlast = powerAllocateXstar;
    end
end
    pstar = xstar;
    Vout = 0;
     VoutUser = zeros(K,1);
    for iUser = 1:K
        Vout = Vout+b2 * pstar' * Mdiag(:,:,1,iUser) * pstar + 1.5 * b4 * norm(pstar' * Mdiag(:,:,1,iUser) * pstar)^2;
        VoutUser(iUser) = b2 * pstar' * Mdiag(:,:,1,iUser) * pstar + 1.5 * b4 * norm(pstar' * Mdiag(:,:,1,iUser) * pstar)^2;
        for iSubbandNumber = 1:subbandNumber
           % sprecoder((iSubbandNumber-1)*Mt+1:(iSubbandNumber-1)*Mt+Mt,1) = pstar( iSubbandNumber) * GainFre(:,iSubbandNumber,1).'/norm(GainFre(:,iSubbandNumber,1));
            sprecoder = xstar;
        if iSubbandNumber ~=1
            Vout = Vout + 3 * b4 * norm(pstar' * Mdiag(:,:,iSubbandNumber,iUser) * pstar)^2;
            VoutUser(iUser) =  VoutUser(iUser) + 3 * b4 * norm(pstar' * Mdiag(:,:,iSubbandNumber,iUser) * pstar)^2;
        end
        end
    end
   
end