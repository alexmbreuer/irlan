%%=======================================================================
%% Copyright 2015 United States Government as represented by Army
%% Research Laboratory.  No copyright is claimed in the United States
%% under Title 17, U.S.Code. All Other Rights Reserved.
%% Authors: Alexander Breuer
%%
%% Distributed under the NASA Open Source Agreement v. 1.3; see
%% accompanying file LICENSE.txt
%%=======================================================================
classdef ThickRestarter < handle
   properties
      dim, wanted
   end

   methods
      function self = ThickRestarter( wanted, dim )
         self.wanted = wanted;
         self.dim = dim;
      end

      function roots = get( self, theta, nconv, varargin )
         roots = theta(1:self.dim - max(self.wanted,floor((3*self.dim + 2*nconv)/5)));
      end
   end
end
