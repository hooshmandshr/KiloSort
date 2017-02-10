function rez = setTemplates(rez)

nt0             = rez.ops.nt0;
rez.ops.nt0min  = ceil(20 * nt0/61);

ops = rez.ops;

rng('default');
rng(1);

Nbatch      = rez.temp.Nbatch;
Nbatch_buff = rez.temp.Nbatch_buff;

Nfilt 	= ops.Nfilt;

ntbuff  = ops.ntbuff;
NT  	= ops.NT;

Nrank   = ops.Nrank;
Th 		= ops.Th;
maxFR 	= ops.maxFR;

Nchan 	= ops.Nchan;

dWU = zeros(nt0, Nchan, Nfilt);
%left window size
lws = ceil(nt0/2)-1;
%right window size
rws = floor(nt0/2);
middle = ceil(size(rez.init_templates, 1)/2);
dWU = single(rez.init_templates(middle-lws:middle+rws, :, :));
% for k = 1:size(dWU, 3)
%     dWU(:, :, k) = (fliplr(dWU(:, :, k)'))';
% end
rez.dWU = dWU;


[W, U, mu, UtU, nu] = decompose_dWU(ops, dWU, Nrank, rez.ops.kcoords);
W0 = W;
W0(NT, 1) = 0;
fW = fft(W0, [], 1);
fW = conj(fW);
rez.fW = fW;


% cid = rez.init_clusterid;
% spt = rez.init_spiketime;
% nspikes = zeros(Nfilt, Nbatch);
% for ibatch = 1:Nbatch
%     from = 1 + (ibatch-1)*NT;
%     to = ibatch*NT;
%     
%     idx = (spt >= from) & (spt < to); 
%     for k = 1:max(cid)
%         nspikes(k, ibatch) = sum(cid(idx)==k);
%     end
% end
% rez.nspikes = nspikes;

end

