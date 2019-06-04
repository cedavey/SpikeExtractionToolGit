% [tlim, vlim] = getTimeAndVoltageLimits(tseries, type)
function [tlim, vlim] = getTimeAndVoltageLimits(tseries, type)
   gettlim = false; getvlim = false;
   if nargin>1 && strcmpi(type,'tlim'), gettlim = true; end
   if nargin>1 && strcmpi(type,'vlim'), getvlim = true; end
   if nargin<2, gettlim = true; getvlim = true; end
   
   if iscell( tseries.data )
      % sometimes some families are in a cell with event based info, while
      % other families are in an array, which means you have to extract
      % mins separately (if 1 family is nuts big but the other isn't)
      if getvlim
         nd   = length( tseries.data );
         dmin = zeros( nd, 1 ); dmax = zeros( size(dmin) );
         for di=1:nd
            % if spikes are made event based the data is kept in a cell
            % array of structs with fields 'spikes' & 'stimes'
            if iscell( tseries.data{di} )
               dmin(di) = min( cellfun(@(d) ternaryOp(isempty(d(:)), 0, min(d.spikes(:))), tseries.data{di} ) );
               dmax(di) = max( cellfun(@(d) ternaryOp(isempty(d(:)), 0, max(d.spikes(:))), tseries.data{di} ) );
            else
               dmin(di) = min( tseries.data{di}(:) );
               dmax(di) = max( tseries.data{di}(:) );
            end
         end
         vlim(1) = min( dmin ) * 1.2;
         vlim(2) = max( dmax ) * 1.2;
      end
      if gettlim
         if ~isempty(tseries.time)
            tlim = [ tseries.time(1) tseries.time(end) ];
         else
            tlim = [tseries.dt tseries.dt*2];
         end
      end
   else
      if getvlim
         vlim = [ nanmin(tseries.data(:)) nanmax(tseries.data(:)) ];
      end
      if gettlim
         tlim = [ tseries.time(1) tseries.time(end) ];
      end
   end
   varargout = cell(nargout, 1);
   if getvlim && gettlim
      varargout{1} = tlim; 
      varargout{2} = vlim; 
   elseif getvlim
      varargout{1} = vlim;
   elseif gettlim
      varargout{1} = tlim;
   end
end
