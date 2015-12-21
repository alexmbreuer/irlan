%%=======================================================================
%% Copyright 2015 United States Government as represented by Army
%% Research Laboratory.  No copyright is claimed in the United States
%% under Title 17, U.S.Code. All Other Rights Reserved.
%% Authors: Alexander Breuer
%%
%% Distributed under the NASA Open Source Agreement v. 1.3; see
%% accompanying file LICENSE.txt
%%=======================================================================
function self = irlan( A, k, x0, dims, tol, maxiter, restarter )
   function advance
      beta = norm(self.r);
      self.Q(:,self.k+1) = self.r/beta;
      self.r = self.A*self.Q(:,self.k+1);
      self.niter = self.niter + 1;
     
      v = self.Q(:,1:self.k+1)'*self.r;
      self.r -= self.Q(:,1:self.k+1)*v;
      v2 = self.Q(:,1:self.k+1)'*self.r;
      self.r -= self.Q(:,1:self.k+1)*v2;
      self.T(self.k+1,self.k) = beta;
      self.T(self.k,self.k+1) = beta;
      self.T(self.k+1,self.k+1) = v(end,end) + v2(end,end);
      self.k = self.k + 1;
   end

   function apply_shift( shifts )
      Tm = self.T;
      S = eye(size(self.T,1));
      for i=1:length(shifts)
         for ii=0:length(self.T)-1
            Tm(end-ii,1:self.k-2-ii) = Tm(1:self.k-2-ii,end-ii) = 0;
	 end
         [W,R] = qr(Tm(1:self.k,1:self.k) - eye(self.k)*shifts(i));
         Tm(1:self.k,1:self.k) = W'*(Tm(1:self.k,1:self.k)*W);
	 S = S*W;
      end

      Q(:,1:self.k) = self.Q(:,1:self.k)*S;

      self.Q = Q;
      self.T = Tm;
      l = length(shifts)-1;
      self.r = self.Q(:,self.k-l)*self.T(self.k-l,self.k-l-1) + self.r*S(self.k,self.k-l-1);
      self.k -= length(shifts);					   
   end

   function check()
      X = self.A*self.Q(:,1:self.k) - self.Q(:,1:self.k)*self.T(1:self.k,1:self.k);
      X(:,self.k) -= self.r;
      assert( norm(X(:,1:self.k),2) < 1e-9*size(X,2)*norm(self.T,2) );
   end

   self.start_time = cputime;
   self.A = A;
   self.Q = zeros(size(A,1),dims);
   self.T = zeros(dims,dims);
   self.mind = k;
   self.niter = 1;
   self.theta_hist = [];
   self.theta_idx = [];
	
   % first iteration
   self.Q(:,1) = x0/norm(x0);
   self.r = self.A*self.Q(:,1);
   self.T(1,1) = self.r'*self.Q(:,1);
   self.r -= self.Q(:,1)*self.T(1,1);
   self.k = 1;
   self.n = dims;

   while( self.k < self.n ) 
      advance();
   end

   [U,Theta] = eig(self.T(1:self.k,1:self.k));
   theta = diag(Theta);
   [theta,perm] = sort(theta);
   U = U(:,perm);
   residuals = abs(U(end,:)*norm(self.r));
   res = sqrt(sum(power( residuals(size(U,1)-self.mind+1:end), 2 )))
   nconv = sum( residuals(size(U,1)-self.mind+1:end) < tol);
   while( nconv < self.mind && self.niter < maxiter )
      shifts = restarter.get( theta, nconv, residuals );
      if( all(shifts == theta(1:size(shifts,1))) )
         self.theta_hist = [self.theta_hist; theta];
         self.theta_idx = [self.theta_idx; size(theta,1)];         
      else
         self.theta_hist = [self.theta_hist; shifts];
         self.theta_idx = [self.theta_idx; -size(shifts,1)];
      end
      stepsize = self.k - k;
      j = 1;
      while j < length(shifts)
         if j + stepsize > length(shifts)
	    apply_shift( shifts(j:end) );
   	    j = length(shifts) + 1;
	 else
            apply_shift( shifts(j:j+stepsize) );
   	    j = j + stepsize;
	 end
         while( self.k < self.n )
     	    advance();
	 end
      end


      [U,Theta] = eig(self.T(1:self.k,1:self.k));
      theta = diag(Theta);
      [theta,perm] = sort(theta);
      U = U(:,perm);
      residuals = abs(U(end,:)*norm(self.r));
      res = sqrt(sum(power( residuals(size(U,1)-self.mind+1:end), 2 )))
      nconv = sum(abs(U(end,size(U,1)-self.mind+1:end)*norm(self.r)) < tol)
   end

   Lambda = theta(size(theta,1)-self.mind+1:end);
   U = self.Q*U(:,size(U,2)-self.mind+1:end);
   self.theta = theta(size(self.Q,2)-self.mind+1:end);
   self.U = U;
   self.end_time = cputime;
end
