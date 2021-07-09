function [corr,asymetry] = correlate_leg_signals(s1_l,s1_r,s2_l,s2_r)
%correlate_leg_signals compute correlation of time signals on two legs

l_corr = mean(normxcorr2(real(s1_l),real(s2_l)));
r_corr = mean(normxcorr2(real(s1_r),real(s2_r)));
corr = mean([l_corr, r_corr]);
asymetry = abs((l_corr-r_corr)/...
                                (l_corr+r_corr));
end

