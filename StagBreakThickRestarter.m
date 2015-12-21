%%=======================================================================
%% Copyright 2015 United States Government as represented by Army
%% Research Laboratory.  No copyright is claimed in the United States
%% under Title 17, U.S.Code. All Other Rights Reserved.
%% Authors: Alexander Breuer
%%
%% Distributed under the NASA Open Source Agreement v. 1.3; see
%% accompanying file LICENSE.txt
%%=======================================================================
classdef StagBreakThickRestarter < handle
   properties
      wanted, dim, w, tau, history, hist_size, hist_valid, hist_idx, min_interval, trcount, thetamin, cheb_degree
   end

   methods
      function self = StagBreakThickRestarter( wanted, dim, w, tau, min_interval, cheb_degree )
         self.wanted = wanted;
	 self.dim = dim;
	 self.w = w;
	 self.tau = tau;
	 self.history = zeros(dim - wanted,w);
	 self.hist_size = size(self.history,1);
	 self.hist_valid = logical(0);
	 self.hist_idx = 1;
	 self.min_interval = min_interval;
	 self.trcount = 0;
	 self.thetamin = [];
	 self.cheb_degree = cheb_degree;
      end

      function roots = get( self, theta, nconv, residuals, varargin )
         if self.hist_valid
	    theta_norm = norm(theta(1:self.hist_size));
	    max_cos = max((self.history'*theta(1:self.hist_size))/theta_norm);
	 else
	    theta_norm = norm(theta(1:self.hist_size));
	    max_cos = -1;
	 end

	 self.history(:,self.hist_idx) = theta(1:self.hist_size)/theta_norm;
	 self.hist_idx = mod(self.hist_idx+1,self.w) + 1;
	 if ~self.hist_valid && self.hist_idx == 1
	    self.hist_valid = logical(1);
	 end

	 if length(self.thetamin) == 0 || self.thetamin(1) > theta(1)
	    self.thetamin(1) = theta(1);
	    self.thetamin(2) = residuals(1);
	 end

	 if 1 - max_cos < self.tau && self.trcount > 2
	    breaking = 1
	    roots = self.cheb_nodes( self.thetamin(1) - self.thetamin(2), self.thetamin(1) );
	    self.trcount = 0;
	 else
            roots = theta(1:self.dim - max(self.wanted,floor((3*self.dim + 2*nconv)/5)));
	    self.trcount = self.trcount + 1;
	 end
      end

      function roots = cheb_nodes( self, a, b )
         roots = zeros(self.cheb_degree,1);
	 for i=1:self.cheb_degree
	    roots(i) = ((1 + cos((pi/2.)*(2*i - 1)/8.))/2.)*(b-a) + a;
	 end
      end
   end
end
