%%=======================================================================
%% Copyright 2015 United States Government as represented by Army
%% Research Laboratory.  No copyright is claimed in the United States
%% under Title 17, U.S.Code. All Other Rights Reserved.
%% Authors: Alexander Breuer
%%
%% Distributed under the NASA Open Source Agreement v. 1.3; see
%% accompanying file LICENSE.txt
%%=======================================================================
classdef Leja < handle
 properties
   maxroots, roots, ab, candidate_points, zprod
 end
 methods
   function self = Leja( maxroots )
      self.maxroots = maxroots;
      self.roots = [];
      self.ab = logical(0);
      self.candidate_points = logical(0);
      self.zprod = [];
   end

   function new_roots = get( self, degree, a, b )
      new_roots = [];

      if ~self.ab
         self.ab = [a b];
	 self.roots = [self.roots a];
	 new_roots = [new_roots a];
	 self.roots = [self.roots b];
	 new_roots = [new_roots b];
	 new_roots = [new_roots a + (b - a)/2.];
	 self.update( a + (b - a)/2. );
      else
         if a < self.ab(1)
	    self.ab(1) = a;
	 end
	 if b > self.ab(2)
	    self.ab(2) = b;
	 end
      end
	 
    for i=1:degree-length(new_roots)
       s = self.candidate_points;

 	 dlo = self.roots(1) - self.ab(1);
 	 clo = self.roots(1) - dlo/2.;
 	 zlo = prod(abs(self.roots - clo));

 	 dhi = self.ab(2) - self.roots(end);
 	 chi = self.roots(end) + dhi/2.;
 	 zhi = prod(abs(self.roots - chi));

 	 if zhi > max(self.zprod)
 	    new_roots = [new_roots chi];
 	    self.update( chi );
 	 elseif zlo > max(self.zprod)
 	    new_roots = [new_roots clo];	    
 	    self.update( clo );
 	 else
 	    [m,idx] = max(self.zprod);
	    
 	    self.candidate_points = [self.candidate_points(1:idx-1) self.candidate_points(idx+1:end)];
 	    self.zprod = [self.zprod(1:idx-1) self.zprod(idx+1:end)];
 	    new_roots = [new_roots s(idx)];
 	    self.update( s(idx) );
       end
    end

    if length(self.roots) > self.maxroots
       self.roots = [self.ab(1) self.ab(2)];
       self.candidate_points = [self.ab(1) + (self.ab(2) - self.ab(1))/2.];
       self.zprod = [prod(abs(self.roots - self.candidate_points))];
    end
   end	 

   function self = update( self, r )
      lidx = lookup( self.roots, r );

      if length(self.zprod) > 0
         self.zprod = self.zprod.*abs( self.candidate_points - r);
      end

      if lidx == 0
         wright = (self.roots(lidx+1) - r)/2.;

   	 self.candidate_points = [self.candidate_points r + wright];
   	 
   	 self.roots = sort([self.roots r]);
   	 self.zprod = [self.zprod prod(abs(self.candidate_points(end) - self.roots))];
      else if lidx == length(self.roots)
         wleft = (r - self.roots(lidx))/2.;

   	 self.candidate_points = [self.candidate_points r - wleft];	          
   	 self.roots = sort([self.roots r]);
   	 self.zprod = [self.zprod prod(abs(self.candidate_points(end) - self.roots))];
      else
   	wleft = (r - self.roots(lidx))/2.;
   	wright = (self.roots(lidx+1) - r)/2.;

   	if self.candidate_points
           self.candidate_points = [self.candidate_points, r - wleft,r + wright];
   	else
           self.candidate_points = [r - wleft,r + wright];
   	end

   	self.roots = sort([self.roots, r]);
   	if length(self.zprod) > 0
   	   len = length(self.candidate_points);
   	   self.zprod = [self.zprod prod(abs(self.candidate_points(len-1) - self.roots)) prod(abs(self.candidate_points(len) - self.roots))];
   	else
   	   len = length(self.candidate_points);
   	   self.zprod = [prod(abs(self.candidate_points(len-1) - self.roots)) prod(abs(self.candidate_points(len) - self.roots))];
   	end	
      end
     end
   end
 end
end
