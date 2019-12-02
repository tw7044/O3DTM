function dtm_values = convert_jd_to_dtm(jd_values)
% CONVERT_JD_TO_DTM allows easy conversion of Julian Date dates to
% Gregorian Date datetimes.
%
% dtm_values = CONVERT_JD_TO_DTM(jd_values)

dtm_values = datetime(jd_values,'convertfrom','juliandate');

end