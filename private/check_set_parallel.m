%function par_state = check_set_parallel(cfg)

v = version('-release');
if str2double(v(1:4)) >= 2014
    v = 'newer';
if str2double(v(1:4)) < 2013
    v = 'older';    
end


if isfield(cfg,'TEparallel') && ft_hastoolbox('DCT');
    if isfield(cfg.TEparallel,'parON') && strcmp(cfg.TEparallel.parON,'yes')
        
        par_state = 1;
        
        switch v
            case {'newer'}
                parallelConfig = parcluster(parallel.defaultClusterProfile);
            case {'2013a' '2013b' 'older'}
                parallelConfig = findResource('scheduler','configuration',defaultParallelConfig);
        end
        max_workers = parallelConfig.ClusterSize;
        if ~isfield(cfg.TEparallel,'workers')
            cfg.TEparallel.workers = max_workers;
        end
        
        switch v
            case {'2013b' 'newer'}
                if  parpool('size')== 0
                    
                    if cfg.TEparallel.workers <= max_workers
                        parpool(cfg.TEparallel.workers)
                    else
                        parpool(max_workers)
                    end
                    
                else if parpool('size') > cfg.TEparallel.workers
                        parpool close;
                        parpool(cfg.TEparallel.workers)
                    end
                end
                
            case {'2013a' 'older'}
                if  matlabpool('size') == 0
                    
                    if cfg.TEparallel.workers <= max_workers
                        matlabpool(cfg.TEparallel.workers)
                    else
                        matlabpool(max_workers)
                    end
                    
                else if matlabpool('size') > cfg.TEparallel.workers
                        matlabpool close;
                        matlabpool(cfg.TEparallel.workers)
                    end
                end
            otherwise
                error('Unknown MATLAB version')
        end
    else
        par_state = 0;
    end
else
    par_state = 0;
end
else
    par_state = 0;
end