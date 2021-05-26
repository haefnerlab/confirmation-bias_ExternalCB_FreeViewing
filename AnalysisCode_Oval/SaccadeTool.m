classdef SaccadeTool < handle
	
	methods ( Access = private )
		function obj = SaccadeTool()
		end
	end

	methods ( Static )
		function sac = Saccade()
			sac.latency		= [];
			sac.duration	= [];
			sac.angle		= [];	% in degrees
			sac.amplitude	= [];
			sac.termiPoints	= single([]);	% 2-by-2: 1st column for start point; 2nd for end point
			sac.velocity	= single([]);
			sac.speed		= single([]);
			sac.peakSpeed	= single([]);
		end

		function [ saccades ] = GetSacs( eyeTrace, samRate, PLOT_FIGURE )
			%% [ saccades ] = GetSacs( eyeTrace, samRate )
			%	analysis the given eye trace and return all saccades detected.
			%	eyeTrace: 2-row matrix in which the first row referring to horizontal eye position
			%			  and the second row designating vertical eye position.
			%	samRate:  sampling rate of eyeTrace.
			%	saccades: a structure array containing properties of all saccades detected.

			if nargin == 1
				samRate = 1000;
				PLOT_FIGURE = 0;
			elseif nargin == 2
				PLOT_FIGURE = 0;
			elseif nargin ~= 3
				disp( 'saccades = GetSacs( eyeTrace[, samRate, PLOT_FIGURE ] )' );
				return;
			end

			if( isempty(eyeTrace) || size(eyeTrace,1) ~= 2 )
				saccades = [];
				disp('eyeTrace must be a 2-row matrix!');
				return;
			end

			maxSamRate = 1000;		% eyeTrace will be resampled with the sampling rate maxSamRate
										% when samRate larger than maxSamRate

			minDur = 0 * MK_CONSTANTS.TIME_UNIT;	% min duration of a saccade
			nDots = size(eyeTrace,2);				% number of sample dots
			if nDots / samRate < minDur				% too short sample duration
				return;
			end
			
			%% resample eyeTrace with the sampling rate maxSamRate
			if samRate > maxSamRate
				step = samRate / maxSamRate;
				nDots = fix( nDots / step );
				for i = 1 : nDots
					eyeTrace(:,i) = mean( eyeTrace( :, fix( (i-1)*step ) + 1 : fix( i*step ) ), 2 );
				end
				eyeTrace( :, nDots+1 : size(eyeTrace,2) ) = [];
				samRate = maxSamRate;
			end
			if PLOT_FIGURE
				figure;
				subplot(221); hold on;
				%plot( eyeTrace' );
			end

			%% deal with high frequency noise with iterated convolution of averaging
			convStep = ceil( max( 0.011 * samRate, 1 ) );
			convFunctor = ones(1,convStep)./convStep;
			eyex = eyeTrace(1,:);
			eyey = eyeTrace(2,:);
			for i=1:100
				%% this "overdone" averaging is equivalent to using a Gaussian functor
				% General model Gauss1:
				% 	f(x) =  a1*exp(-((x-b1)/c1)^2)
				% 	Coefficients (with 95% confidence bounds):
				% 		a1 =     0.01261  (0.01261, 0.01262)	( \sigma	= 31.64  )
				% 		b1 =         501  (501, 501)			( \mu		= 501 )
				% 		c1 =        44.7  (44.69, 44.7)
				%
				% 	Goodness of fit:
				% 		SSE: 2.288e-08
				% 		R-square: 1
				% 		Adjusted R-square: 1
				% 		RMSE: 4.788e-06

				eyexl = [ ones( 1, fix(convStep/2) ) * eyex(1), eyex, ones( 1, fix(convStep/2) ) * eyex(nDots) ];
				eyexl = conv( eyexl, convFunctor, 'same' );
				eyex = eyexl( fix(convStep/2) + 1: nDots + fix(convStep/2) );
				eyeyl = [ ones( 1, fix(convStep/2) ) * eyey(1), eyey, ones( 1, fix(convStep/2) ) * eyey(nDots) ];
				eyeyl = conv( eyeyl, convFunctor, 'same' );
				eyey = eyeyl( fix(convStep/2) + 1: nDots + fix(convStep/2) );
				if i == 1
					eyeTrace = [ eyex; eyey ];
				end
			end
			smoothEyeTrace = [ eyex; eyey ];
			if PLOT_FIGURE
				plot( eyeTrace' );
				plot( smoothEyeTrace' );
			end

			%% calculate velocity
			velocity = gradient( smoothEyeTrace, 1/samRate );	% 1/samRate specifies the space between points
			convStep = ceil( max( 0.005 * samRate, 1 ) );
			convFunctor = ones(1,convStep)./convStep;
			for i=1:10
				velocity(1,:) = conv( velocity(1,:), convFunctor, 'same' ); % 'same' to get the central part
				velocity(2,:) = conv( velocity(2,:), convFunctor, 'same' ); % 'same' to get the central part
			end
			[ polVel.angle, polVel.speed ] =  cart2pol( velocity(1,:), velocity(2,:) );	% velocity in polar coordinate

			%% calculate acceleration
			acceleration = gradient( velocity, 1/samRate );
			aclFromSpeed = gradient( polVel.speed, 1/samRate );
			convStep = ceil( max( 0.005 * samRate, 1 ) );
			convFunctor = ones(1,convStep)./convStep;
			for i=1:10
				acceleration(1,:) = conv( acceleration(1,:), convFunctor, 'same' );
				acceleration(2,:) = conv( acceleration(2,:), convFunctor, 'same' );
			end
			[ polAccel.angle, polAccel.acceleration ] = cart2pol( acceleration(1,:), acceleration(2,:) );
			
			if PLOT_FIGURE
				subplot(222); hold on;
				plotyy( 1:nDots, velocity(1,:), 1:nDots, acceleration(1,:) );
				subplot(223);
				plotyy( 1:nDots, velocity(2,:), 1:nDots, acceleration(2,:) );
				subplot(224);
				plotyy(  1:nDots, polVel.speed, 1:nDots, aclFromSpeed );
			end


			%% look for saccades
			% pick out saccade candidates using a speed threshold
			minSpeed = 2.1;		% min speed for a saccade
			mark = zeros(1,nDots);
			mark( find( polVel.speed > minSpeed ) ) = 1;
			mark(1) = 0;
			mark(end) = 0;

			iSpeedTroughs = find( aclFromSpeed(1:end-1) <= 0 & aclFromSpeed(2:end) > 0 );	% trough position of polVel.Speed
			mark( iSpeedTroughs ) = 0;	% for a small minSpeed, two adjacent saccades are very easy to be recognized as one;
			mark( iSpeedTroughs+1 ) = 0;	% the trough position is used to seperate two adjacent saccades in such a case

			if PLOT_FIGURE
				subplot(221);
				plot( mark );
			end

			mark = mark(2:end) - mark(1:end-1);
			iBounds = [ find( mark == 1 ); find( mark == -1 ) ];	% first row for start positions and second for end ones	
			if isempty(iBounds)
				saccades = [];
				return;
			end
			if iBounds( 1, 1 ) < 0.005 * samRate	% too close to trial start
				iBounds( :, 1 ) = [];
		    end
		    if isempty(iBounds)
				saccades = [];
				return;
		    end
			if ( nDots - iBounds( 2, end) ) / samRate < 0.005	% too close to trial end
				iBounds( :, end ) = [];
		    end
		    if isempty(iBounds)
				saccades = [];
				return;
		    end

			if PLOT_FIGURE
				plot(iBounds(1,:),ones(1,size(iBounds,2))*2,'r.');
				plot(iBounds(2,:),ones(1,size(iBounds,2))*2,'r.');
			end

			% peakShift measures the ratio of the distance between the peak position(speed) and start position
			%to that between the end position and the peak position; for a saccade, it's usually close to 1
			iSpeedPeaks = zeros(1,size(iBounds,2));
			iLefts = zeros(1,size(iBounds,2));
			iRights = zeros(1,size(iBounds,2));
			for i = 1:length(iSpeedPeaks)
				%index = find( aclFromSpeed( iBounds(1,i):iBounds(2,i)-1 ) >= 0 & aclFromSpeed( iBounds(1,i)+1:iBounds(2,i) ) < 0 );
				index = find( polVel.speed( iBounds(1,i):iBounds(2,i) ) == max( polVel.speed( iBounds(1,i):iBounds(2,i) ) ) );
				iSpeedPeaks(i) = iBounds(1,i) + index( ceil(length(index)/2) ) - 1;
				index = find( polVel.speed( iBounds(1,i):iSpeedPeaks(i) ) <= 1/3 * polVel.speed(iSpeedPeaks(i)), 1, 'last' );
				if ~isempty(index)
					iLefts(i) = iBounds(1,i) + index;
				else
					iLefts(i) = iBounds(1,i);
				end
				%iLefts(i) = iBounds(1,i);
				index = find( polVel.speed( iSpeedPeaks(i):iBounds(2,i) ) <= 1/3 * polVel.speed(iSpeedPeaks(i)), 1, 'first' );
				if ~isempty(index)
					iRights(i) = iSpeedPeaks(i) + index - 2;
				else
					iRights(i) = iBounds(2,i);
				end
				%iRights(i) = iBounds(2,i);
			end
			if PLOT_FIGURE
				plot(iSpeedPeaks(1,:),ones(1,length(iSpeedPeaks))*3,'r.');
			end
			peakShift = ( iSpeedPeaks - iLefts ) ./ ( iRights - iSpeedPeaks );
			peakShift(  peakShift == 1/0  ) = 0;
			peakShiftThresh = 0.25;
			iBounds( : , abs( peakShift - 1 ) > peakShiftThresh ) = [];
			iSpeedPeaks( abs( peakShift - 1 ) > peakShiftThresh ) = [];
			peakShift( abs( peakShift - 1 ) > peakShiftThresh ) = [];
			if PLOT_FIGURE
				plot(iBounds(1,:),ones(1,size(iBounds,2))*2.2,'g.');
				plot(iBounds(2,:),ones(1,size(iBounds,2))*2.2,'g.');
				plot(iSpeedPeaks(1,:),ones(1,length(iSpeedPeaks))*2.8,'g.');
			end

			if isempty(iBounds)
				saccades = [];
				return;
			end

			%% refine start positions and end positions	
			velocity = gradient( eyeTrace, 1/samRate );	% 1/samRate specifies the space between points
			convStep = ceil( max( 0.005 * samRate, 1 ) );
			convFunctor = ones(1,convStep)./convStep;
			for i = 1 : 5
				velocity(1,:) = conv( velocity(1,:), convFunctor, 'same' ); % 'same' to get the central part
				velocity(2,:) = conv( velocity(2,:), convFunctor, 'same' ); % 'same' to get the central part
			end
			[ polVel.angle, polVel.speed ] =  cart2pol( velocity(1,:), velocity(2,:) );	% velocity in polar coordinate
			
			minSpeedThresh = 3;
			lowl2peak = zeros(1,size(iBounds,2));
			lowr2peak = zeros(1,size(iBounds,2));
			saccades( size(iBounds,2) ) = SaccadeTool.Saccade();
			for i = 1 : size( iBounds, 2 )
				minSpeed = max( minSpeedThresh, 0.05 * max( polVel.speed( iBounds(1,i) : iBounds(2,i) ) ) );	% 0.05 of the peak speed		
				index = find( polVel.speed( iBounds(1,i) : iSpeedPeaks(i) ) <= minSpeed, 1, 'last' );
				if ~isempty(index)
					iBounds(1,i) = iBounds(1,i) + index;
				end
				index = find( polVel.speed( iSpeedPeaks(i) : iBounds(2,i) ) <= minSpeed, 1, 'first' );
				if ~isempty(index)
					iBounds(2,i) = iSpeedPeaks(i) - 2 + index;
				end

				index = find( polVel.speed( iBounds(1,i):iBounds(2,i) ) == max( polVel.speed( iBounds(1,i):iBounds(2,i) ) ) );
				if isempty(index)
					saccades(i).duration = 0;
					continue;
				end

				lowl2peak(i) = polVel.speed( iBounds(1,i) ) / max( polVel.speed( iBounds(1,i) : iBounds(2,i) ) );
				lowr2peak(i) = polVel.speed( iBounds(2,i) ) / max( polVel.speed( iBounds(1,i) : iBounds(2,i) ) );

				iSpeedPeaks(i) = iBounds(1,i) + index( ceil(length(index)/2) ) - 1;

				saccades(i).latency = ( iBounds(1,i) - 1 ) / samRate * MK_CONSTANTS.TIME_UNIT;
				saccades(i).velocity = velocity( : , iBounds(1,i) : iBounds(2,i) );
				saccades(i).speed = polVel.speed( iBounds(1,i) : iBounds(2,i) );
				saccades(i).peakSpeed = max( saccades(i).speed );
				saccades(i).duration = ( iBounds(2,i) - iBounds(1,i) ) / samRate * MK_CONSTANTS.TIME_UNIT;
				[ saccades(i).angle, saccades(i).amplitude ] = cart2pol( eyeTrace(1,iBounds(2,i)) - eyeTrace(1,iBounds(1,i)), eyeTrace(2,iBounds(2,i)) - eyeTrace(2,iBounds(1,i)) );
				saccades(i).angle = saccades(i).angle * 180 / pi;
				saccades(i).termiPoints = [ eyeTrace(1,iBounds(1,i)), eyeTrace(1,iBounds(2,i)); eyeTrace(2,iBounds(1,i)), eyeTrace(2,iBounds(2,i)) ];
			end
			iBounds( :, [saccades.duration] <= 0 ) = [];
			iSpeedPeaks( [saccades.duration] <= 0 ) = [];
			lowl2peak( [saccades.duration] <=0 ) = [];
			lowr2peak( [saccades.duration] <=0 ) = [];
			saccades( [saccades.duration] <= 0 ) = [];
			% peakShift = ( iSpeedPeaks - iBounds(1,:) ) ./ ( iBounds(2,:) - iSpeedPeaks )
			% size(peakShift)
			% size(iSpeedPeaks)
			% peakShiftThresh = 0.45;
			% saccades( abs( peakShift - 1 ) > peakShiftThresh ) = [];

			%% delete saccades with high ratio of low speed to peak speed ( very likely to be influenced by low frequency noise )
			saccades( lowl2peak >= 0.24 | lowr2peak >= 0.24 ) = [];		% this should be equivalent to use a threshold of 12.5 (velocity threshold 3) for the peak velocity

			%% delete saccades with short intervals
			min_interval = 0.02 * MK_CONSTANTS.TIME_UNIT;
			index = ( [ saccades(2:end).latency ] - ( [ saccades(1:end-1).latency ] + [ saccades(1:end-1).duration ] ) ) < min_interval;
			index = [ index, 0 ] | [ 0, index ];
			saccades( index ) = [];

			%% delete saccades with improper durations
			min_duration = 0.01 * MK_CONSTANTS.TIME_UNIT;
			max_duration = 0.15 * MK_CONSTANTS.TIME_UNIT;
			saccades( [ saccades.duration ] >= max_duration ) = [];
			saccades( [ saccades.duration ] <= min_duration ) = [];

			if PLOT_FIGURE
				plot(iBounds(1,:),ones(1,size(iBounds,2))*2.2,'m.');
				plot(iBounds(2,:),ones(1,size(iBounds,2))*2.2,'m.');
				plot(polVel.speed./50);
				%plot(polVel.speed);
			end

			%figure
			%hist(polVel.speed,min(polVel.speed):0.2:max(polVel.speed))
		end

		function map = GetCalibrationMap( isForward, set_points, recorded_points )
			%% isForward:		true means map from recorded points to set points, vice versa
			%  set_points:		2-column/row array specifies set points
			%  recorded_points:	2-column/row array specifies recorded points
			%  map:				specifies the relationship between a set points and a recorded points
			map = [];
			if( isempty(set_points) ||...
				size( size(set_points), 2 ) ~= size( size(recorded_points), 2 ) ||...
				any( size(set_points) ~= size(recorded_points) ) )
				return;
			end
			if( size(set_points,2) > 2 )
				set_points = set_points';
				recorded_points = recorded_points';
			end

			if( isForward )
				input_points = double( set_points );
				base_points = double( recorded_points );
			else
				input_points = double( recorded_points );
				base_points = double( set_points );
			end

			switch size(set_points,1)
				case 1
					return;
				case 2
					map = cp2tform( input_points, base_points, 'nonreflective similarity' );
				case 3
					map = cp2tform( input_points, base_points, 'affine' );
				case { 4, 5 }
					map = cp2tform( input_points, base_points, 'projective' );
				case { 6, 7, 8, 9 }
					map = cp2tform( input_points, base_points, 'polynomial', 2 );
				case { 10, 11, 12, 13, 14 }
					map = cp2tform( input_points, base_points, 'polynomial', 3 );
				otherwise
					map = cp2tform( input_points, base_points, 'polynomial', 4 );
			end


		end

		function [ x, y ] = Calibrate( map, isForward, x, y )
			%% isForward: true means map from recorded points to set points, vice versa
			if( isempty(map) ) return; end
			if(isForward)	[x,y] = tforminv( map, double(x), double(y) );
			else			[x,y] = tformfwd( map, double(x), double(y) );
			end
			x = single(x);
			y = single(y);
		end
	end

end