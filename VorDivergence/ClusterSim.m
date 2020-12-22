classdef ClusterSim < handle
    %genClusters Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N = 1000; % number of particles
        NClusters = 10; % Number of clusters
        MeanNPerCluster = 200; % average points per cluster
        ClusterSigma = 2; % STD of the points from their expected positions
        SetType = 'Diffusive'; % determines which Cluster Style is Simulated
        ROI = [128 128];
        RandomSeed = 8675309;
        UseSeed = false;
        NoiseDensity = 0.25; % per sim unit^2
    end
    
    % Gaussian only properties
    properties
        
    end
    
     % Diffusion only prpoerties
    properties
        DiffCte = 0.5; % Diffusion Constante
        TrackLength = 100; % number of samples to take
    end
    
    % Tubule only properties
    properties(Transient = true)
        Tubule1 = [0,0,128,128];
        Tubule2 = [0,128,128,0];
        ClusterLongSigma = 2;
        ClusterShortSigma = 0.5;
    end
    
    properties(SetAccess=protected)
        Types = {'Gaussian','Diffusive','Tubule'}
        Clusters;
        CoordinateList;
        TypeGenerated;
        ClusterSpawnPoints
    end
    
    methods
        function obj = ClusterSim(varargin)
            %ClusterSim Constructor
            if obj.UseSeed
               rng(obj.RandomSeed) 
            end
        end
        
        function obj = simulateData(obj)
            switch obj.SetType
                case 'Gaussian'
                    obj.getGaussianSim;
                case 'Diffusive'
                    obj.getDiffusiveSim;
                case 'Tubule'
                    obj.getTubuleSim;
            end
        end
        
        function obj = getGaussianSim(obj)
            % generate cluster spawn points
            ClusterCenters = rand(obj.NClusters,2)*7/8; % add a buffer of 1/8 ROI
            ClusterCenters = [ClusterCenters(:,1)*obj.ROI(1) ClusterCenters(:,2)*obj.ROI(2)];
            % cluster points
            PperCluster = poissrnd(obj.MeanNPerCluster,[obj.NClusters,1]);
            RandomVectors = obj.ClusterSigma*randn(sum(PperCluster),2);
            GaussCluster = cell(obj.NClusters,1);
            indexStart = 1;
            for ii = 1:obj.NClusters
               GaussCluster{ii} = ClusterCenters(ii,:) ...
                   + RandomVectors(indexStart:indexStart+PperCluster(ii)-1,:);
               indexStart = indexStart+PperCluster(ii);
            end
            % background noise
            NoisePositions = obj.getBackgroundNoise;
            % output data
            obj.Clusters = GaussCluster;
            obj.CoordinateList = [cell2mat(GaussCluster); NoisePositions];
            obj.TypeGenerated = 'Gaussian';
            obj.ClusterSpawnPoints = num2cell(ClusterCenters,2);
        end
        
        function obj = getDiffusiveSim(obj)
            % generate cluster spawn points
            ClusterStarts = rand(obj.NClusters,2)*2/3; % add a buffer of 1/3 ROI
            ClusterStarts = [ClusterStarts(:,1)*obj.ROI(1) ClusterStarts(:,2)*obj.ROI(2)];
            tracks = cell(obj.NClusters,1);
            for ii = 1:obj.NClusters
                tracks{ii} = cumsum(sqrt(2*obj.DiffCte)*randn(obj.TrackLength,2),1);
                tracks{ii} = [tracks{ii}(:,1)+ClusterStarts(ii,1),tracks{ii}(:,2)+ClusterStarts(ii,2)];                
            end
            % cluster points
            PperCluster = poissrnd(obj.MeanNPerCluster,[obj.NClusters,1]);
            RandomTrackPos = randi(obj.TrackLength,sum(PperCluster),1);
            RandomVectors = obj.ClusterSigma*randn(sum(PperCluster),2);
            DiffuseCluster = cell(obj.NClusters,1);
            indexStart = 1;
            for ii = 1:obj.NClusters
               DiffuseCluster{ii} = tracks{ii}(RandomTrackPos...
                   (indexStart:indexStart+PperCluster(ii)-1),:)...
                   +RandomVectors(indexStart:indexStart+PperCluster(ii)-1,:);
               indexStart = indexStart+PperCluster(ii);
            end
            % background noise
            NoisePositions = obj.getBackgroundNoise;
            % output data
            obj.Clusters = DiffuseCluster;
            obj.CoordinateList = [cell2mat(DiffuseCluster); NoisePositions];
            obj.TypeGenerated = 'Diffusive';
            obj.ClusterSpawnPoints = tracks;
        end
        
        function obj = getTubuleSim(obj)
            % generate cluster spawn points
            Slope1 = (obj.Tubule1(2)-obj.Tubule1(4))/(obj.Tubule1(1)-obj.Tubule1(3));
            Slope2 = (obj.Tubule2(2)-obj.Tubule2(4))/(obj.Tubule2(1)-obj.Tubule2(3));
            Intercept1 = obj.Tubule1(2)-obj.Tubule1(1)*Slope1;
            Intercept2 = obj.Tubule2(2)-obj.Tubule2(1)*Slope2;
            AngleTubule1 = atan2((obj.Tubule1(4)-obj.Tubule1(2)),(obj.Tubule1(3)-obj.Tubule1(1)));
            AngleTubule2 = atan2((obj.Tubule2(4)-obj.Tubule2(2)),(obj.Tubule2(3)-obj.Tubule2(1)));
            halfPoint = round(obj.NClusters/2);
            CenterCandidates = rand(obj.NClusters,1);
            %TauCenters = zeros(obj.NClusters,2);
            xPos = abs(obj.Tubule1(3)-obj.Tubule1(1))*CenterCandidates(1:halfPoint)+obj.Tubule1(1);
            xPos = [xPos; abs(obj.Tubule2(3)-obj.Tubule2(1))*CenterCandidates(halfPoint+1:end)+obj.Tubule2(1)];
            yPos = [xPos(1:halfPoint)*Slope1+Intercept1; xPos(halfPoint+1:end)*Slope2+Intercept2];
            TauCenters = [xPos yPos];
            % cluster points
            PperCluster = poissrnd(obj.MeanNPerCluster,[obj.NClusters,1]);
            RandomVectors = [obj.ClusterShortSigma*randn(sum(PperCluster),1),...
                obj.ClusterLongSigma*randn(sum(PperCluster),1)];
            % rotate random vectors to align with tubule positions
            Tubule1Count = sum(PperCluster(1:halfPoint));
            firstTubuleX = RandomVectors(1:Tubule1Count,1)*cos(AngleTubule1)...
                + RandomVectors(1:Tubule1Count,2)*sin(AngleTubule1);
            firstTubuleY = -RandomVectors(1:Tubule1Count,1)*sin(AngleTubule1)...
                + RandomVectors(1:Tubule1Count,2)*cos(AngleTubule1);
            secondTubuleX = RandomVectors(Tubule1Count+1:end,1)*cos(AngleTubule2)...
                + RandomVectors(Tubule1Count+1:end,2)*sin(AngleTubule2);
            secondTubuleY = -RandomVectors(Tubule1Count+1:end,1)*sin(AngleTubule2)...
                + RandomVectors(Tubule1Count+1:end,2)*cos(AngleTubule2);
            TubuleDeviations = [firstTubuleX firstTubuleY; secondTubuleX secondTubuleY];
            % assign to clusters
            TubuleCluster = cell(obj.NClusters,1);
            indexStart = 1;
            for ii = 1:obj.NClusters
               TubuleCluster{ii} = TauCenters(ii,:) ...
                   + TubuleDeviations(indexStart:indexStart+PperCluster(ii)-1,:);
               indexStart = indexStart+PperCluster(ii);
            end
            % background noise
            NoisePositions = obj.getBackgroundNoise;
            % output data
            obj.Clusters = TubuleCluster;
            obj.CoordinateList = [cell2mat(TubuleCluster); NoisePositions];
            obj.TypeGenerated = 'Tubule';
            obj.ClusterSpawnPoints = num2cell(TauCenters,2);
        end
        
        function NoisePositions = getBackgroundNoise(obj)
            expectedNoiseCount = obj.NoiseDensity*prod(obj.ROI);
            NoiseCount = poissrnd(expectedNoiseCount);
            NoisePositions = [rand(NoiseCount,1)*obj.ROI(1) rand(NoiseCount,1)*obj.ROI(2)];
        end
    end
end

