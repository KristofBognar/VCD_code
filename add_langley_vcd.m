function [avg_vcd] = add_langley_vcd(rcd_S, avg_vcd)

N = size(avg_vcd.day);
for i = 1:1:N(1)
   TF = (avg_vcd.day(i) == rcd_S.day) & (avg_vcd.ampm(i) == rcd_S.ampm);
   avg_vcd.langley_vcd(i,:) = rcd_S.m(TF,:);
   avg_vcd.langley_vcd_err(i,:) = rcd_S.m_1sigma(TF,:);
end
