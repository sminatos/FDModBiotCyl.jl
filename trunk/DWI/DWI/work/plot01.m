clear all
close all

filename='test.8ht'
A=load(filename);
ns=A(8);
nrec=length(A)/(ns+8)
B=reshape(A,ns+8,nrec);
data=B(9:end,:);
%figure;imagesc(data)



% $$$     //   header --> struct TMP_HEAD @include/ovsp.h
% $$$     //                int   seq;
% $$$     //                int   rec_dist;      // (cm)
% $$$     //                int   shot_dist;     // (cm)
% $$$     //                int   gain;          // (db)
% $$$     //                int   rec_component; // V(1), H1(2), H2(3)
% $$$     //                float delay;       // (microsec)
% $$$     //                float sampling;    // (microsec)
% $$$     //                int   data_length;
% $$$     //         
% $$$     //  wavefrom_header[i].data_length = nt / n_skip;

dt=A(7)*1E-6;
tvec=[0:dt:(ns-1)*dt];
figure;plot(tvec,-data(:,1))
%tmpa=axis;axis([0 3E-3 tmpa(3) tmpa(4)])
figure;plot(tvec,-data(:,1)/max(abs(data(:,1))))

