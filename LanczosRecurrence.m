%%=======================================================================
%% Copyright 2015 United States Government as represented by Army
%% Research Laboratory.  No copyright is claimed in the United States
%% under Title 17, U.S.Code. All Other Rights Reserved.
%% Authors: Alexander Breuer
%%
%% Distributed under the NASA Open Source Agreement v. 1.3; see
%% accompanying file LICENSE.txt
%%=======================================================================
classdef LanczosRecurrence
  properties
     A, Q, T, r, n, mind, k
  end
  methods
     function self = LanczosRecurrence( A, x0, n, mind )
        self.A = A;
     	self.Q = zeros(size(A,1),n);
	self.T = zeros(n,n);
	self.mind = mind;
	
	% first iteration
	self.Q(:,1) = x0/norm(x0);
	self.r = self.A*self.Q(:,1);
	self.T(1,1) = self.r'*self.Q(:,1);
	self.r -= self.Q(:,1)*self.T(1,1);
	self.k = 1;
	self.n = 11;
     end

     function self = advance( self ) 
        beta = norm(self.r);
	self.T(self.k+1,self.k) = self.T(self.k,self.k+1) = beta;
        self.Q(:,self.k+1) = self.r/beta;
	self.r = self.A*self.Q(:,self.k+1);
	
	v = self.Q(:,1:self.k+1)'*self.r;
	self.r -= self.Q(:,1:self.k+1)*v;
	v2 = self.Q(:,1:self.k+1)'*self.r;
	self.r -= self.Q(:,1:self.k+1)*v2;
	self.T(self.k+1,self.k+1) = v(end,end) + v2(end,end);
	self.k = self.k + 1;
     end

     function self = apply_shift( self, s )
        [W,R] = qr(self.T(1:self.k,1:self.k) - eye(self.k)*s);
	self.T(1:self.k,1:self.k) = W'*self.T(1:self.k,1:self.k)*W;
	self.Q(:,1:self.k) = self.Q(:,1:self.k)*W;

	self.r = self.Q(:,k)*self.T(self.k,self.k-1) - self.r*W(self.k,self.k-1);
	self.k -= 1;
     end
  end
end
