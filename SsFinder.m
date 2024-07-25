classdef SsFinder
    methods(Static)
        function best_match = findPss(samples,max_freq_shift,samples_per_symb)
            best_match=struct("max_abs",-1);
            for id=0:2
                pss=[PssGenerator.generatePssByCellInfo(id) zeros(1,samples_per_symb-127)];
                for freq=(0:max_freq_shift)+56
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
            % figure
            % sss=[SssGenerator.generateSssByCellInfo(NId1+NId2) zeros(1,samples_per_symb-127)];
            % corr=xcorr(area,ifft(circshift(sss,56+kssb)));
            % plot(abs(corr));
        end
    end
end