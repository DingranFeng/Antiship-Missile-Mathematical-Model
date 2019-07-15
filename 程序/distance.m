function d=distance(lon1,lat1,lon2,lat2)
%东经1，北纬1，东经2，北纬2 (度)
RE=6300; %地球半径 km
d=RE*acos(cos(pi/180*lat1)*cos(pi/180*lat2)*cos(pi/180*(lon1-lon2))+sin(pi/180*lat1)*sin(pi/180*lat2)); %km
end

