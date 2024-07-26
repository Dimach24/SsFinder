classdef SsFinder
    methods(Static)
        function best_match = findPss(samples,from_freq,to_freq,samples_per_symb)
            best_match=struct("max_abs",-1);
            for id=0:2
                pss=[PssGenerator.generatePssByCellInfo(id) zeros(1,samples_per_symb-127)];
                for freq=(from_freq:to_freq)+56
                    corr=xcorr(samples,ifft(circshift(pss,freq)));
                    if best_match.max_abs<max(abs(corr))
                        best_match.NId2=id;
                        best_match.freq=freq;
                        best_match.corr=corr;
                        best_match.max_abs=max(abs(corr));
                    end
                end
            end
            [~,best_match.lags]=xcorr(samples,pss);
            best_match.kSSB=best_match.freq-56;
        end
        function [NId1,max_corr] = checkSss(samples, pss_index, kssb, NId2, samples_per_symb)
            area=samples(pss_index+2*samples_per_symb:pss_index+3*samples_per_symb);
            max_corr=-1;
            for id=(0:335)*3+NId2
                sss=[SssGenerator.generateSssByCellInfo(id) zeros(1,samples_per_symb-127)];
                corr=xcorr(area,ifft(circshift(sss,56+kssb)));
                current=max(abs(corr));
                if current>max_corr
                    max_corr=current;
                    NId1=floor(id/3);
                end
            end
        end
        function peaks_i=findPeaks(samples,level)
            arguments
                samples
                level = 100
            end
            trs=(level<=samples);
            peaks_i=[];
            max_val=-1;
            max_i=-1;
            is_seq=false;
            for i=1:length(samples)
                if trs(i)
                    if ~is_seq
                        is_seq=true;
                        max_val=samples(i);
                        max_i=i;
                    else
                        if (samples(i))>max_val
                            max_i=i;
                            max_val=samples(i);
                        end
                    end
                else
                    if is_seq
                        is_seq=false;
                        peaks_i(end+1)=max_i;
                    end
                end
            end
        end
        function [NCellId,kSSB,t_index,cropped]=processSignalByPeakNo(...
                samples,...
                from_freq,...
                to_freq,...
                samples_per_symb,...
                peak_No)
            match=SsFinder.findPss(samples,from_freq,to_freq,samples_per_symb);
            peaks=match.lags(SsFinder.findPeaks(abs(match.corr),110));
            NId1=SsFinder.checkSss(samples,peaks(1),match.kSSB,match.NId2,samples_per_symb);
            NCellId=NId1*3+match.NId2;
            kSSB=match.kSSB;
            t_index=peaks(peak_No);
            cropped=samples(t_index:end);
        end
    end
end