a = load('freqsortn.txt');
nlist = a(:,1);
llist = a(:,2);
flist = a(:,3);
num = length(nlist);
fmax = max(flist);
lmax = max(llist);
l = [];
f = [];
n = 0;
for i = 1:num
    if nlist(i) == n
        l(end+1) = llist(i);
        f(end+1) = flist(i);
    else
        n = n+1;
        plot(l,f,'k.-','LineWidth',0.01,'MarkerSize',4);
        hold on
        l = [];
        f = [];
    end
end
ylim([0 fmax+5])
xlim([0 lmax+5])
xlabel('angular degree $ l $','interpreter','latex');
ylabel('eigenfrequency (mHz)')
hold off
