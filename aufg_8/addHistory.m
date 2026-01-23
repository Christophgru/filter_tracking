function Hist = addHistory(Hist, val)

global HIST_SIZE;
if(isempty(Hist))
  for i=1:HIST_SIZE
     Hist(:,i) = val;
  end
end
Hist = circshift(Hist',1)';       
Hist(:,1) = val;

end

