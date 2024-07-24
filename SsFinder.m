classdef SsFinder
    methods(Static)
        function best_match = findPss(samples,max_freq_shift)
            best_match=struct("max_abs",-1);
            for id=0:2
                pss=[PssGenerator.generatePssByCellInfo(id) zeros(1,1920-127)];
                for freq=(0:max_freq_shift)+56
                    corr=xcorr(samples,ifft(circshift(pss,freq)));
                    if best_match.max_abs<max(abs(corr))
                        best_match.id=id;
                        best_match.freq=freq;
                        best_match.corr=corr;
                        best_match.max_abs=max(abs(corr));
                    end
                end
            end
            [~,best_match.lags]=xcorr(samples,pss);
            best_match.freq=best_match.freq-56;
        end
    end
end