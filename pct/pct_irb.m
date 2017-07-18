function pct_irb(V, meta, irb_file)
%PCT_IRB converts CTP data into IRB format
%
%   Ruogu Fang 07/06/2013
%   Advanced Multimedia Processing (AMP) Lab, Cornell University
%
%   USAGE:  PCT_IRB(V, META, IRB_FILE);
%
%   PRE:
%       V           - CTP input data [Y x X x Z x T]
%       META        - A meta file information [1x9] vector
%       IRB_FILE    - Filename for data in IRB format
%

slice = meta(1);
AIFy = meta(2);
AIFx = meta(3);
VOFy = meta(4);
VOFx = meta(5);
PRE = meta(6);
POST = meta(7);

V = permute(double(squeeze(V(:,:,slice,:))),[3 1 2]);

save(irb_file,'V','AIFy','AIFx','VOFy','VOFx','PRE','POST');

end





    
    