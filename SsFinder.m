classdef SsFinder
    methods(Static)
        function [sample_indexes, k_SSB, NId1] = findPss(samples,max_freq_shift,peak_finder)
            corr_res=zeros(3,max_freq_shift+1,length(samples));
            for id=0:2
                pss=[PssGenerator.generatePssByCellInfo(id) zeros(max_freq_shift)];
                for freq=0:max_freq_shift
                    corr_res(id+1,freq+1,:)=xcorr(ifft(circshift(pss,freq)),samples);
                end
            end
            absolute=abs(corr_res);
            m=zeros(3);
            f_max=zeros(3);
            for id=1:3
                [m(id),f_max(id)]=max(absolute(id));
                m(id)=max(m(id));
            end
            [~,id]=max(m);
            array_with_peaks=absolute(id,f_max(id));
            NId1=id-1;
            k_SSB=f_max-57; % -56 -1
            sample_indexes=peak_finder(array_with_peaks);
            
        end
    end
end