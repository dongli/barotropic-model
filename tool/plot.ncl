load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    start_time = 0
    end_time = 2
    time_step = 1
    file_prefix = "output"
    file_root = "../build/Debug"

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; i = i+"+time_step+" )); do "+ \
                    "    printf '"+file_root+"/"+file_prefix+".%3.3d.nc\n' $i; done")

    wks = gsn_open_wks("pdf", "plot")

    res_common = True
    res_common@gsnFrame = False
    res_common@gsnDraw = False
    res_common@gsnSpreadColors = True
    res_common@cnFillOn = True
    res_common@cnLinesOn = False
    res_common@cnLevelSelectionMode = "ManualLevels"

    res1 = res_common
    res1@cnMinLevelValF = 80000
    res1@cnMaxLevelValF = 90000
    res1@cnLevelSpacingF = 500

    res2 = res_common

    res3 = res_common

    res_panel = True

    do i = 0, dimsizes(fs)-1
        system("echo time step: "+i)
        f = addfile(fs(i), "r")

        plot1 = gsn_csm_contour_map(wks, f->gh(0,:,:), res1)
        plot2 = gsn_csm_contour_map(wks, f->u(0,:,:), res2)
        plot3 = gsn_csm_contour_map(wks, f->v(0,:,:), res3)

        gsn_panel(wks, (/plot1,plot2,plot3/), (/3,1/), res_panel)
    end do

end