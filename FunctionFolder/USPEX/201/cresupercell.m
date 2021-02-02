function [newlattice,newcoord,newtypes,neworders,numIons]=cresupercell(lattice, coordinates, atomtype, order , cell1, cell2)
global ORG_STRUC      
%%This program is to covert the cell1 to the cell2
%%The basic idea is like this, let take the example of coversion from 2*2 to 3*3
%%firstly creat supercell 6*6 based on 2*2, and then cut the supercell as 4
%%regions, finally take one region randomly.
       newlattice=lattice;
       newcoord=[];
       neworders=[];
       newtypes=[];
       numIons=[];
   if size(coordinates,1)==0   
       return
   else 
%%%%%1:create supercell

       super(1)=lcm(cell1(1),cell2(1));
       super(2)=lcm(cell1(2),cell2(2));
       num1=super(1)/cell1(1);
       num2=super(2)/cell1(2);
       item=1;

       for i=1:size(coordinates,1)
          for j=0:num1-1
	       for k=0:num2-1                
                   coord(item,1)=coordinates(i,1)/num1;
                   coord(item,2)=coordinates(i,2)/num2;
                   coord(item,3)=coordinates(i,3);

                   coord(item,1)=coord(item,1)+j/num1;
                   coord(item,2)=coord(item,2)+k/num2;
                   S_type(item)=atomtype(i);                  
                   S_order(item)=order(i);       
                   item=item+1;
	         	end
            end
           
       end
       Slattice(1,:)=lattice(1,:)*num1;
       Slattice(2,:)=lattice(2,:)*num2;
       Slattice(3,:)=lattice(3,:);

%%%%%%2: chose the region
      num3=super(1)/cell2(1);
      num4=super(2)/cell2(2);
      coord=coord-floor(coord);
%%%%%%%%%%%%%%how many atoms in each region
max1=0;
pick1=[];
pick2=[];
     for i=1:num3
         for j=1:num4
             ioncount=0;
             for k=1:item-1
                 if ( ((i-1)/num3-0.1) < coord(k,1) )&&  ( coord(k,1) < (i/num3+0.1) )
                     if(  ((j-1)/num4-0.1) < coord(k,2))  && ( coord(k,2) < (j/num4+0.1) )
                         ioncount=ioncount+1;
                     end
                 end
             end
             if ioncount>max1
                max1=ioncount;
                pick1=i;
                pick2=j;
             end      
         end
     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      num_Ions=0;
      newlattice(1,:)=Slattice(1,:)/num3;
      newlattice(2,:)=Slattice(2,:)/num4;
      newlattice(3,:)=Slattice(3,:);
      for i=1:item-1
          if( ((pick1-1)/num3-0.1) < coord(i,1) ) && ( coord(i,1) < (pick1/num3+0.1) )
              if(  ((pick2-1)/num4-0.1) < coord(i,2))  && ( coord(i,2) < (pick2/num4+0.1) )
                  num_Ions=num_Ions+1;
                  newcoord(num_Ions,:)=coord(i,:)*Slattice/newlattice;
		  newtypes(num_Ions)=S_type(i);
                  neworders(num_Ions)=S_order(i);
              end
	  end
      end
      for i=1:size(ORG_STRUC.atomType,2)
          if(isempty(newtypes))
              numIons=[];
          else
              numIons(i)=sum(size(find(newtypes==ORG_STRUC.atomType(i)),2));
          end
      end
   end
