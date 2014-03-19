load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

    start_time = 0
    end_time = 23940
    time_step = 180
    file_prefix = "output"
    file_root = "../build/Debug"
    only_plot_gh = True

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; i = i+"+time_step+" )); do "+ \
                    "    printf '"+file_root+"/"+file_prefix+".%5.5d.nc\n' $i; done")

    wks = gsn_open_wks("pdf", "plot")

    res_common = True
    res_common@gsnMaximize = True
    res_common@gsnFrame = False
    res_common@gsnDraw = False
    res_common@mpCenterLonF = 180
    res_common@gsnPolar = "NH"
    res_common@gsnSpreadColors = True
    res_common@cnLineLabelsOn = False
    res_common@cnLineThicknessF = 3.0
    res_common@cnLevelSelectionMode = "ManualLevels"

    res1 = res_common
    res1@cnLevelSelectionMode = "ManualLevels"
    res1@cnMinLevelValF = 0
    res1@cnMaxLevelValF = 50
    res1@cnLevelSpacingF = 5

    res2 = res_common
    res2@cnLevelSelectionMode = "ManualLevels"
    res2@cnMinLevelValF = -40
    res2@cnMaxLevelValF = 40
    res2@cnLevelSpacingF = 5

    res3 = res_common
    res3@cnMinLevelValF = 80000
    res3@cnMaxLevelValF = 100000
    res3@cnLevelSpacingF = 1000

    res_panel = True

    do i = 0, dimsizes(fs)-1
        system("echo time step "+i)
        f = addfile(fs(i), "r")

        if (only_plot_gh) then
            res3@gsnRightString = "Day "+cd_string(f->time, "%j")
            plot3 = gsn_csm_contour_map(wks, f->gh(0,:,:), res3)
            draw(plot3)
            frame(wks)
        else
            plot1 = gsn_csm_contour_map(wks, f->u(0,:,:), res1)
            plot2 = gsn_csm_contour_map(wks, f->v(0,:,:), res2)
            plot3 = gsn_csm_contour_map(wks, f->gh(0,:,:), res3)
            gsn_panel(wks, (/plot1,plot2,plot3/), (/2,2/), res_panel)
        end if
    end do

end