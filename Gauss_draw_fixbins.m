function [histo_gap,histo_condu] = Gauss_draw_fixbins(Conductance, bins)

[m,~] = size(Conductance);
histo_gap = linspace(0,2,bins);
histo_condu = zeros(1,bins);
for i=1:bins-1
    for j=1:m
        if (Conductance(j,1) >  histo_gap(i))&&(Conductance(j,1)<histo_gap(i+1))
            histo_condu(1,i)= histo_condu(1,i)+1;
        end
    end
end