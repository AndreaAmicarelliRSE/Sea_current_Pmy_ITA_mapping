program Sea_current_Pmy_ITA_mapping
!-------------------------------------------------------------------------------
! "Sea current Pmy ITA mapping v.1.0" 
! Copyright 2016 (RSE SpA)
! "Sea current Pmy ITA mapping v.1.0" authors and email contact are provided on 
! the documentation file.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Description: This tool reads the ".nc" files on marine current yearly average     
!              power from Sea_current_Py_mapping and provides:
!                 1) the multiple year average of the specific power flow;
!                 2) the 2D multiple year average of the specific power flow 
!                    (highest values over the depth);
!                 3) same as 2), but for the coastal regions;
!                 4) same as 3), but for Italian coasts.
!              The output is written in binary ".dat", ".vtk" and ".txt" file 
!              formats. 
!-------------------------------------------------------------------------------  
!------------------------
! Modules
!------------------------ 
use netcdf
!------------------------
! Declarations
!------------------------
implicit none
! The 3D yearly matrices are 72 x 677 x 253 matrices 
integer,parameter :: ndepth=72,nlat=253,nlon=677
! Italian approximate boundaries: lon_ITA_min=7., lat_ITA_min=36.
integer,parameter :: ilon_ITA_min=208,ilat_ITA_min=87,nlon_ITA=210,nlat_ITA=167
integer :: i,j,k,IOstatus,k1,k2,ndata,nfile
real :: P_max,P_min,distance
character(100) FILE_NAME
real :: lon(nlon),lat(nlat),depth(ndepth)
real :: power_y(nlon,nlat,ndepth),power_my(nlon,nlat,ndepth)
real :: power_zstar(nlon,nlat),power_coast(nlon,nlat)
real :: power_ITA(nlon_ITA,nlat_ITA),power_coast_den(nlon,nlat)
integer :: n_power_my(nlon,nlat,ndepth)
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
P_max = -1000.
P_min = 1000.
write(*,*)                                                                     &
"This tool reads the ".nc" files on marine current yearly average ",           &
"power from Sea_current_Py_mapping and provides: ",                            &
"   1) the multiple year average of the specific power flow; ",                &
"   2) the 2D multiple year average of the specific power flow (highest ",     &
"      values over the depth); ",                                              &
"   3) same as 2), but for the coastal regions; ",                             &
"   4) same as 3), but for Italian coasts. ",                                  &
"The output is written in binary ".dat", ".vtk" and ".txt" file formats. "
do j=1,nlat
   do i=1,nlon
      power_zstar(i,j) = -999.
      power_coast(i,j) = -999.
      power_coast_den(i,j) = 0.
      do k=1,ndepth
         power_my(i,j,k) = 0.
         n_power_my(i,j,k) = 0
      enddo
   enddo
enddo
nfile = 0
!------------------------
! Statements
!------------------------
! Start input reading and pre-processing
write(*,*) "   Start input reading and pre-processing"
! Reading coordinates
open(1,file='lon.dat',form='unformatted')
read(1) lon
close (1)
open(1,file='lat.dat',form='unformatted')
read (1) lat
close (1)
open(1,file='depth.dat',form='unformatted')
read(1) depth
close (1)
! Reading file name list
write(*,*) "   Start reading file name list "
open(1,file='file_name_list.txt')
! Cycle over files
do
   read(1,'(a)',IOSTAT=IOstatus) FILE_NAME
   if (IOstatus/=0) exit
! Reading yearly average files
   write(*,*) "   Reading file ",FILE_NAME
   nfile = nfile + 1
   open(2,file=FILE_NAME,form='unformatted')
   read(2) power_y
   close(2)
! Contribution to the mutiple year average P; maxima and minima yearly values
   do k=1,ndepth    
      do j=1,nlat
         do i=1,nlon
            if (power_y(i,j,k)/=-999.) then
               power_my(i,j,k) = power_my(i,j,k) + power_y(i,j,k)
               n_power_my(i,j,k) = n_power_my(i,j,k) + 1
               if (power_y(i,j,k)>P_max) P_max = power_y(i,j,k)
               if (power_y(i,j,k)<P_min) P_min = power_y(i,j,k)
            endif
         enddo
      enddo
   enddo
enddo
close(1)
write(*,*) "   End reading file name list "
write(*,*) "   P_max(W/m^2) ",P_max," P_min(W/m^2) ",P_min
write(*,*) "   End reading and pre-processing yearly average data"
! End input reading and pre-processing
! Computation of the multiple year average of the specific power flow 
write(*,*) "2) Start computation of the multiple year P "
do k=1,ndepth    
   do j=1,nlat
      do i=1,nlon
         if (n_power_my(i,j,k)>=nfile/2) then
            power_my(i,j,k) = power_my(i,j,k) / n_power_my(i,j,k)
            else
               power_my(i,j,k) = -999.
         endif
      enddo
   enddo
enddo
write(*,*) "   End computation of the multiple year P "
! Computation of the 2D power field (maxima multiple year average over the 
! depth)
write(*,*) "   Start Computation of the 2D power field (maxima of 2) over ",   &
           "      the depth) "
do k=1,ndepth    
   do j=1,nlat
      do i=1,nlon
         if (power_my(i,j,k)>power_zstar(i,j)) then
            power_zstar(i,j) = power_my(i,j,k)
         endif
      enddo
   enddo
enddo
write(*,*) "   End Computation of the 2D power field (maxima of 2) over ",     &
           "      the depth) "
!Test start
! 4) Extraction of the coastal region power values (cells close to mainland; P=-1: missing value for marine non-coastal regions)
!  print *,"4) Start Extraction of the coastal region power values "
!  do i = 2, (nlon-1)
!     do j = 2, (nlat-1)
!          if ((power_zstar(i+1,j)==-999.).or.(power_zstar(i-1,j)==-999.).or.(power_zstar(i,j+1)==-999.).or.(power_zstar(i,j-1)==-999.).or. &
!              (power_zstar(i+1,j+1)==-999.).or.(power_zstar(i+1,j-1)==-999.).or.(power_zstar(i-1,j+1)==-999.).or.(power_zstar(i-1,j-1)==-999.)) then
!             power_coast(i,j)=power_zstar(i,j)
!             else
!                if (power_zstar(i,j)/=-999.) power_coast(i,j)=-1.
!          endif
!     end do
!  end do
!  print *,"   End Extraction of the coastal region power values "
! 4) Extrapolation into mainland to fill the gaps with GIS (cells close to mainland; P=-1: missing value for marine non-coastal regions)
!Test end
write(*,*) "   Start Extrapolation into mainland to fill the gaps with GIS "
do j=4,(nlat-3)
   do i=4,(nlon-3)
      if (power_zstar(i,j)==-999.) then
         do k1=1,7
            do k2=1,7
               if (power_zstar(i-4+k1,j-4+k2)/=-999.) then
                  distance = sqrt(((-4 + k1) * 5.2) ** 2 + ((-4 + k2) * 6.9)   &
                             ** 2)
                  if (power_coast(i,j)==-999.) power_coast(i,j) = 0.
                  power_coast(i,j) = power_coast(i,j) +                        &
                                     power_zstar(i-4+k1,j-4+k2) /              &
                                     (distance ** 5)
                  power_coast_den(i,j) = power_coast_den(i,j) + 1. /           &
                                         (distance ** 5) 
               endif
            enddo
         enddo
         if (power_coast_den(i,j)/=0.) power_coast(i,j) = power_coast(i,j) /   &
            power_coast_den(i,j)   
         else
            power_coast(i,j) = power_zstar(i,j)
      endif
   enddo
enddo
write(*,*) "   End Extrapolation into mainland to fill the gaps with GIS "
! Extraction of the italian coastal region power values 
write(*,*) "   Start Extraction of the italian coastal region power values "
do j=ilat_ITA_min,nlat
   do i=ilon_ITA_min,(ilon_ITA_min+nlon_ITA-1)
      power_ITA(i-ilon_ITA_min+1,j-ilat_ITA_min+1) = power_coast(i,j)  
   enddo
enddo
write(*,*) "   End extraction of the italian coastal region power values "
! Output writing (binary ".dat", ASCII ".vtk" and ASCII ".txt")
write(*,*) "   Start output writing (binary .dat, ASCII .vtk and ASCII .txt) "
! Write binary ".dat" output
open(2,file='power_my.dat',form='unformatted')
write(2) power_my
close(2)
open(2,file='power_zstar.dat',form='unformatted')
write(2) power_zstar
close(2)
open(2,file='power_coast.dat',form='unformatted')
write(2) power_coast
close(2)
open(2,file='power_ITA.dat',form='unformatted')
write(2) power_ITA
close(2)
! Write ".vtk" output
! Write 3D multiple year average P
open(3,file='power_my.vtk')
write(3,'(a)') '# vtk DataFile Version 3.0'
write(3,'(a)') 'POWER_MY'
write(3,'(a)') 'ASCII'
write(3,'(a)') 'DATASET RECTILINEAR_GRID'
write(3,'(a,3(i5,1X))') 'DIMENSIONS ',nlon,nlat,ndepth
ndata = nlat * nlon * ndepth
write(3,'(a,i10,a)') 'X_COORDINATES ',nlon,' float '
do i=1,nlon
   write(3,'(F12.4)',ADVANCE='NO') lon(i)
enddo
write(3,'(/a,i10,a)') 'Y_COORDINATES ',nlat,' float '
do i=1,nlat
   write(3,'(F12.4)',ADVANCE='NO') lat(i)
enddo
write(3,'(/a,i10,a)') 'Z_COORDINATES ',ndepth,' float '
do i=1,ndepth
    write(3,'(F12.4)',ADVANCE='NO') depth(i)
enddo
write(3,'(/a,i10)') 'POINT_DATA ',ndata
write (3,'(a)') 'SCALARS A float '
write (3,'(a)') 'LOOKUP_TABLE default ' 
do k=1,ndepth
   do j=1,nlat
      do i=1,nlon
         write (3,'(F12.4)') power_my(i,j,k)
      enddo
   enddo
enddo  
close(3) 
! Write 2D P_zstar
open(3,file='power_zstar.vtk')
write(3,'(a)') '# vtk DataFile Version 3.0'
write(3,'(a)') 'POWER_ZSTAR'
write(3,'(a)') 'ASCII'
write(3,'(a)') 'DATASET RECTILINEAR_GRID'
write(3,'(a,2(i5,1X))') 'DIMENSIONS ',nlon,nlat,1
ndata = nlon * nlat
write(3,'(a,i10,a)') 'X_COORDINATES ',nlon,' float '
do i=1,nlon
   write(3,'(F12.4)',ADVANCE='NO') lon(i)
enddo
write(3,'(/a,i10,a)') 'Y_COORDINATES ',nlat,' float '
do i=1,nlat
   write(3,'(F12.4)',ADVANCE='NO') lat(i)
enddo
write(3,'(a)') 'Z_COORDINATES 1 float'
write(3,'(a)') '0. '
write(3,'(/a,i10)') 'POINT_DATA ',ndata
write(3,'(a)') 'SCALARS A float '
write(3,'(a)') 'LOOKUP_TABLE default ' 
do j=1,nlat
   do i=1,nlon
      write(3,'(F12.4)') power_zstar(i,j)
   enddo
enddo
close(3) 
! Write 2D P_coast
open(3,file='power_coast.vtk')
write(3,'(a)') '# vtk DataFile Version 3.0'
write(3,'(a)') 'POWER_COAST'
write(3,'(a)') 'ASCII'
write(3,'(a)') 'DATASET RECTILINEAR_GRID'
write(3,'(a,2(i5,1X))') 'DIMENSIONS ',nlon,nlat,1
write(3,'(a,i10,a)') 'X_COORDINATES ',nlon,' float '
do i=1,nlon
   write(3,'(F12.4)',ADVANCE='NO') lon(i)
enddo
write(3,'(/a,i10,a)') 'Y_COORDINATES ',nlat,' float '
do i=1,nlat
   write(3,'(F12.4)',ADVANCE='NO') lat(i)
enddo
write(3,'(a)') 'Z_COORDINATES 1 float '
write(3,'(a)') '0. '
write(3,'(/a,i10)') 'POINT_DATA ',ndata
write(3,'(a)') 'SCALARS A float '
write(3,'(a)') 'LOOKUP_TABLE default ' 
do j=1,nlat
   do i=1,nlon
      write(3,'(F12.4)') power_coast(i,j)
   enddo
enddo
close(3) 
! Write 2D P_ITA
open(3,file='power_ITA.vtk')
write(3,'(a)') '# vtk DataFile Version 3.0'
write(3,'(a)') 'POWER_ITA'
write(3,'(a)') 'ASCII'
write(3,'(a)') 'DATASET RECTILINEAR_GRID'
write(3,'(a,2(i5,1X))') 'DIMENSIONS ',nlon_ITA,nlat_ITA,1
ndata = nlon_ITA * nlat_ITA
write(3,'(a,i10,a)') 'X_COORDINATES ',nlon_ITA,' float '
do i=1,nlon_ITA
   write(3,'(F12.4)',ADVANCE='NO') lon(i+ilon_ITA_min-1)
enddo
write(3,'(/a,i10,a)') 'Y_COORDINATES ',nlat_ITA,' float '
do i=1,nlat_ITA
   write(3,'(F12.4)',ADVANCE='NO') lat(i+ilat_ITA_min-1)
enddo
write(3,'(a)') 'Z_COORDINATES 1 float '
write(3,'(a)') '0. '
write(3,'(/a,i10)') 'POINT_DATA ',ndata
write(3,'(a)') 'SCALARS A float '
write(3,'(a)') 'LOOKUP_TABLE default ' 
do j=1,nlat_ITA
   do i=1,nlon_ITA
      write(3,'(F12.4)') power_ITA(i,j)
   enddo
enddo
close(3) 
! Write 2D P_ITA in ".txt" file
open(4,file='power_ITA.txt')
write (4,'(a)') '         lon         lat         P(W/m^2) '
do i=1,nlon_ITA
   do j=1,nlat_ITA
      write(4,"(1X,3(F12.4,1X))") lon(i+ilon_ITA_min-1),lat(j+ilat_ITA_min-1), &
         power_ITA(i,j)
   enddo
enddo
close(4)
write(*,*) "   End output writing (binary ".dat" and ASCII ".vtk") "
write(*,*) "*** Sea_current_Pmy_ITA_mapping has terminated "
!------------------------
! Deallocations
!------------------------
end program Sea_current_Pmy_ITA_mapping

