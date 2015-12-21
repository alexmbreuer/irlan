%%=======================================================================
%% Copyright 2015 United States Government as represented by Army
%% Research Laboratory.  No copyright is claimed in the United States
%% under Title 17, U.S.Code. All Other Rights Reserved.
%% Authors: Alexander Breuer
%%
%% Distributed under the NASA Open Source Agreement v. 1.3; see
%% accompanying file LICENSE.txt
%%=======================================================================
classdef LejaRestarter < handle
   properties
      poly, deg, aidx, bidx
   end
   methods
      function self = LejaRestarter( maxroots, deg, aidx, bidx )
         self.poly = Leja( maxroots );
	 self.deg = deg;
	 self.aidx = aidx;
	 self.bidx = bidx;
      end

      function roots = get( self, theta, varargin )
         roots = self.poly.get( self.deg, theta(self.aidx), theta(self.bidx) );
      end
   end
end
