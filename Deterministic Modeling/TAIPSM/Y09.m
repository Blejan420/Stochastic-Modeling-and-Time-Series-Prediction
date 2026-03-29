%
%	Data file Y09.M
%
% Monthly average of weekly working time in industry and services, in the USA, between 1979 and 1989. 
%

	y=[35.3 35.8 35.8 35.8 35.9 36.3 36.4 36.3 36.1 36.1 35.9 36.2 35.1 ...
	   35.3 35.8 35.8 35.7 36.2 36.3 36.2 35.9 35.9 35.8 36.8 35.2 35.4 ...
	   35.7 35.1 35.5 35.9 36.0 36.0 35.8 35.7 35.6 35.9 35.1 35.1 35.2];
	y=[y ...
	   35.0 35.0 35.3 35.3 35.5 35.3 35.3 35.3 35.7 35.1 35.0 35.2 35.2 ...
	   35.2 35.4 35.6 35.6 35.1 35.2 35.1 35.2 33.9 34.8 34.8 34.7 34.8 ...
	   35.0 35.2 35.2 34.8 34.7 34.7 35.0 34.6 34.2 34.7 34.7 34.9 35.2];
	y=[y ...
	   35.4 35.4 35.3 35.3 35.1 35.5 34.9 34.9 35.0 35.2 35.1 35.4 35.5 ...
	   35.4 35.3 35.0 35.0 35.4 34.6 34.5 34.8 34.7 34.9 35.2 35.1 35.2 ...
	   35.0 34.9 34.8 35.2];

	Ts = 1 ; 
	unit = 'Time [months]' ; 
	ntime  = 0:(length(y)-1) ;
	label = 'Monthly average of weekly working time in industry and services, in the USA (1979-1989).' ; 
	yunit = 'Monthly average of weekly working time [hours]' ;

%
% Author:  Dan Stefanoiu
% Date:    30.10.1992 
% Updated: 25.07.2007
%