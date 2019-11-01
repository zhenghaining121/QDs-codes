function [Gauss_peak_ori,sigma] = Gauss_peak_fit(M_ori,bins_ori)

Mx_ori = zeros(bins_ori, bins_ori);
My_ori = zeros(bins_ori, bins_ori);
Mymax_ori = zeros(1, bins_ori);
Gauss_peak_ori = zeros(1, bins_ori);
sigma = zeros(1,bins_ori);

for i = 1:bins_ori
    Mx_ori(:,i)=((1:bins_ori)/bins_ori).*ones(1,bins_ori);
    My_ori(:,i)=fliplr(M_ori(:,i))/fliplr(max(M_ori(:,i)));
    Mymax_ori(i)=max(M_ori(:,i));
    fun=fittype('A*exp(-(x-mu)^2/(2*sigma^2))');
    [cf,~]=fit(Mx_ori(:,i),My_ori(:,i),fun,'Start',[0.5 0.5 0.5]);
    ki_ori=linspace(1,bins_ori,50000);
    xi_ori=interp1(1:bins_ori,Mx_ori(:,i),ki_ori,'nearest');
    sigma(1,i) = cf.sigma;
    Yi_ori=cf.A*exp(-(xi_ori-cf.mu).^2/(2*cf.sigma^2))*Mymax_ori(i);
    Gauss_peak_ori(i)=bins_ori*cf.mu;
end
