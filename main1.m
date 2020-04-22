%�򻯰�ABC
% �˹���Ⱥ�㷨����Ѱ��
% ABC�㷨
clc,clear,close all
warning off
format longG

% ABC ������ʼ��
maxiter = 20;  % ��������
sizepop = 10;  % ��Ⱥ����
popmin1 = -1;  popmax1 = 1; % x1
popmin2 = -1;  popmax2 = 1; % x2
trail(1:sizepop,1) = 0;     % δ�ҵ����Ž�ĵ�������
limit = sizepop;            % ��β������

% ��ʼ����Ⱥ
for i=1:sizepop
    x1 = popmin1+(popmax1-popmin1)*rand;
    x2 = popmin2+(popmax2-popmin2)*rand;
    pop(i,:) = [x1, x2];
    fitness(i) = fun([x1,x2]);  % ��Ӧ�Ⱥ���ֵ--Ŀ�꺯��ֵ--��СĿ�꺯��ֵ
end

% ��¼һ������ֵ
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   %ȫ�����
fitnesszbest=bestfitness; %ȫ�������Ӧ��ֵ

% ����Ѱ��
for i=1:maxiter
   
    % ���۷忪ʼ����
    for j=1:sizepop
        % ѡ����۵ĸ���
        x2y = randi(2);   % ����δ֪��
        % ѡ����������Ⱥ
        neighbour = randi(sizepop);  
        % neighbour != j
        while( neighbour==j  )
            neighbour = randi(sizepop);
        end
        
        % ��Ⱥ���£��⣩
        tempx = pop(j,:);  % ��ǰ�Ľ�
        tempx(x2y) = pop(j,x2y) + (pop(j,x2y)-pop(neighbour,x2y))*(rand-0.5)*2;
%         tempx(x2y) = pop(j,x2y) + (pop(j,x2y)-pop(neighbour,x2y))*rand;
        
        % x1  Խ������
        if tempx(1)>popmax1
            tempx(1)=popmax1;
        end
        if tempx(1)<popmin1
            tempx(1)=popmin1;
        end
        % x2  Խ������
        if tempx(2)>popmax2
            tempx(2)=popmax2;
        end
        if tempx(2)<popmin2
            tempx(2)=popmin2;
        end
        
        % ��Ӧ�ȸ���
        temp_fitness = fun( tempx );  % ��ǰ����Ӧ�Ⱥ���ֵ--Ŀ�꺯��ֵ--��СĿ�꺯��ֵ
        % �Ƚ�  �����Ƚ�
        if temp_fitness<fitness(j)
            fitness(j) = temp_fitness;
            pop(j,:) = tempx;
        end
        if temp_fitness<bestfitness
            bestfitness = temp_fitness;
            zbest =  tempx;
        end
    end
   
    % �۲��
    % �������
    prob = 0.9*fitness./max(fitness) + 0.1;
    for j=1:sizepop
        if(rand<prob(j))
            
            % ѡ����۵ĸ���
            x2y = randi(2);   % ����δ֪��
            % ѡ����������Ⱥ
            neighbour = randi(sizepop);  
            % neighbour != j
            while( neighbour==j  )
                neighbour = randi(sizepop);
            end

            % ��Ⱥ���£��⣩
            tempx = pop(j,:);  % ��ǰ�Ľ�
            tempx(x2y) = pop(j,x2y) + (pop(j,x2y)-pop(neighbour,x2y))*(rand-0.5)*2;
%             tempx(x2y) = pop(j,x2y) + (pop(j,x2y)-pop(neighbour,x2y))*rand;

            % x1  Խ������
            if tempx(1)>popmax1
                tempx(1)=popmax1;
            end
            if tempx(1)<popmin1
                tempx(1)=popmin1;
            end
            % x2  Խ������
            if tempx(2)>popmax2
                tempx(2)=popmax2;
            end
            if tempx(2)<popmin2
                tempx(2)=popmin2;
            end

            % ��Ӧ�ȸ���
            temp_fitness = fun( tempx );  % ��ǰ����Ӧ�Ⱥ���ֵ--Ŀ�꺯��ֵ--��СĿ�꺯��ֵ
            % �Ƚ�  �����Ƚ�
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
   
    % ���忪ʼ����
    [maxTrial, index] = max(trail);
    index = index(1);
    if(maxTrial(1)>limit)
        x1 = popmin1+(popmax1-popmin1)*rand;
        x2 = popmin2+(popmax2-popmin2)*rand;
        pop(index,:) = [x1, x2];
        fitness(index) = fun([x1,x2]);  % ��Ӧ�Ⱥ���ֵ--Ŀ�꺯��ֵ--��СĿ�꺯��ֵ
    end
   
    fitness_iter(i) = bestfitness;
end

disp('���Ž�')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
axis tight