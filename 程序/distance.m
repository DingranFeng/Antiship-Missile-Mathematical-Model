function d=distance(lon1,lat1,lon2,lat2)
%����1����γ1������2����γ2 (��)
RE=6300; %����뾶 km
d=RE*acos(cos(pi/180*lat1)*cos(pi/180*lat2)*cos(pi/180*(lon1-lon2))+sin(pi/180*lat1)*sin(pi/180*lat2)); %km
end

