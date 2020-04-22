%简化版ABC
% 人工蜂群算法函数寻优
% ABC算法
clc,clear,close all
warning off
format longG

% ABC 参数初始化
maxiter = 20;  % 迭代次数
sizepop = 10;  % 种群数量
popmin1 = -1;  popmax1 = 1; % x1
popmin2 = -1;  popmax2 = 1; % x2
trail(1:sizepop,1) = 0;     % 未找到更优解的迭代次数
limit = sizepop;            % 拖尾最大次数

% 初始化种群
for i=1:sizepop
    x1 = popmin1+(popmax1-popmin1)*rand;
    x2 = popmin2+(popmax2-popmin2)*rand;
    pop(i,:) = [x1, x2];
    fitness(i) = fun([x1,x2]);  % 适应度函数值--目标函数值--最小目标函数值
end

% 记录一组最优值
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   %全局最佳
fitnesszbest=bestfitness; %全局最佳适应度值

% 迭代寻优
for i=1:maxiter
   
    % 采蜜峰开始工作
    for j=1:sizepop
        % 选择采蜜的个体
        x2y = randi(2);   % 两个未知数
        % 选择相连的种群
        neighbour = randi(sizepop);  
        % neighbour != j
        while( neighbour==j  )
            neighbour = randi(sizepop);
        end
        
        % 种群更新（解）
        tempx = pop(j,:);  % 当前的解
        tempx(x2y) = pop(j,x2y) + (pop(j,x2y)-pop(neighbour,x2y))*(rand-0.5)*2;
%         tempx(x2y) = pop(j,x2y) + (pop(j,x2y)-pop(neighbour,x2y))*rand;
        
        % x1  越界限制
        if tempx(1)>popmax1
            tempx(1)=popmax1;
        end
        if tempx(1)<popmin1
            tempx(1)=popmin1;
        end
        % x2  越界限制
        if tempx(2)>popmax2
            tempx(2)=popmax2;
        end
        if tempx(2)<popmin2
            tempx(2)=popmin2;
        end
        
        % 适应度更新
        temp_fitness = fun( tempx );  % 当前的适应度函数值--目标函数值--最小目标函数值
        % 比较  个体间比较
        if temp_fitness<fitness(j)
            fitness(j) = temp_fitness;
            pop(j,:) = tempx;
        end
        if temp_fitness<bestfitness
            bestfitness = temp_fitness;
            zbest =  tempx;
        end
    end
   
    % 观察峰
    % 计算概率
    prob = 0.9*fitness./max(fitness) + 0.1;
    for j=1:sizepop
        if(rand<prob(j))
            
            % 选择采蜜的个体
            x2y = randi(2);   % 两个未知数
            % 选择相连的种群
            neighbour = randi(sizepop);  
            % neighbour != j
            while( neighbour==j  )
                neighbour = randi(sizepop);
            end

            % 种群更新（解）
            tempx = pop(j,:);  % 当前的解
            tempx(x2y) = pop(j,x2y) + (pop(j,x2y)-pop(neighbour,x2y))*(rand-0.5)*2;
%             tempx(x2y) = pop(j,x2y) + (pop(j,x2y)-pop(neighbour,x2y))*rand;

            % x1  越界限制
            if tempx(1)>popmax1
                tempx(1)=popmax1;
            end
            if tempx(1)<popmin1
                tempx(1)=popmin1;
            end
            % x2  越界限制
            if tempx(2)>popmax2
                tempx(2)=popmax2;
            end
            if tempx(2)<popmin2
                tempx(2)=popmin2;
            end

            % 适应度更新
            temp_fitness = fun( tempx );  % 当前的适应度函数值--目标函数值--最小目标函数值
            % 比较  个体间比较
            if temp_fitness<fitness(j)
                fitness(j) = temp_fitness;
                pop(j,:) = tempx;
            end
            if temp_fitness<bestfitness
                bestfitness = temp_fitness;
                zbest =  tempx;
            else
                trail(j) = trail(j)+1;
            end
        end
    end
   
    % 侦察峰开始工作
    [maxTrial, index] = max(trail);
    index = index(1);
    if(maxTrial(1)>limit)
        x1 = popmin1+(popmax1-popmin1)*rand;
        x2 = popmin2+(popmax2-popmin2)*rand;
        pop(index,:) = [x1, x2];
        fitness(index) = fun([x1,x2]);  % 适应度函数值--目标函数值--最小目标函数值
    end
   
    fitness_iter(i) = bestfitness;
end

disp('最优解')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
axis tight