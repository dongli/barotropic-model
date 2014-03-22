load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

    start_time = 0
    end_time = 3600
    time_step = 60
    file_prefix = "output"
    file_root = "../build/Debug"
    only_plot_gh = False
    plot_contour_line = True
    plot_wind_vector = True

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; i = i+"+time_step+" )); do "+ \
                    "    printf '"+file_root+"/"+file_prefix+".%5.5d.nc\n' $i; done")

    wks = gsn_open_wks("pdf", "plot")

    gsn_define_colormap(wks, "rainbow")

    res_common = True
    res_common@gsnMaximize = True
    res_common@gsnFrame = False
    res_common@gsnDraw = False
    res_common@mpCenterLonF = 180
    res_common@gsnPolar = "NH"

    if (plot_contour_line) then
        res_contour_line = True
        res_contour_line@cnLineLabelsOn = False
        res_contour_line@cnLineThicknessF = 3.0
        res_contour_line@cnLevelSelectionMode = "ManualLevels"
    else
        res_contour_fill = True
        res_contour_fill@cnFillOn = True
        res_contour_fill@cnLinesOn = False
        res_contour_line@cnLevelSelectionMode = "ManualLevels"
        res_contour_fill@gsnSpreadColors = True
    end if

    if (plot_wind_vector) then
        res1 = True
        copy_VarAtts(res_common, res1)
        res1@gsnLeftString = "velocity"
        res1@vcRefLengthF = 0.01
        res1@vcRefMagnitudeF = 10
    else
        res1 = True
        copy_VarAtts(res_common, res1)
        if (plot_contour_line) then
            copy_VarAtts(res_contour_line, res1)
        else
            copy_VarAtts(res_contour_fill, res1)
        end if
        res1@cnMinLevelValF = -50
        res1@cnMaxLevelValF = 50
        res1@cnLevelSpacingF = 5

        res2 = True
        copy_VarAtts(res_common, res2)
        if (plot_contour_line) then
            copy_VarAtts(res_contour_line, res2)
        else
            copy_VarAtts(res_contour_fill, res2)
        end if
        res2@cnMinLevelValF = -50
        res2@cnMaxLevelValF = 50
        res2@cnLevelSpacingF = 5
    end if

    res3 = True
    copy_VarAtts(res_common, res3)
    if (plot_contour_line) then
        copy_VarAtts(res_contour_line, res3)
    else
        copy_VarAtts(res_contour_fill, res3)
    end if
    res3@gsnLeftString = "geopotential height"
    res3@cnMinLevelValF = 70000
    res3@cnMaxLevelValF = 100000
    res3@cnLevelSpacingF = 1000

    res_panel = True

    do i = 0, dimsizes(fs)-1
        system("echo plot file "+fs(i))
        f = addfile(fs(i), "r")
        gh = f->gd+f->ghs
        copy_VarCoords(f->gd, gh)

        if (only_plot_gh) then
            res3@gsnRightString = "Day "+cd_string(f->time, "%j")
            plot3 = gsn_csm_contour_map(wks, gh(0,:,:), res3)
            draw(plot3)
            frame(wks)
        else
            if (plot_wind_vector) then
                plot1 = gsn_csm_vector_map(wks, f->u(0,:,:), f->v(0,:,:), res1)
                plot3 = gsn_csm_contour_map(wks, gh(0,:,:), res3)
                gsn_panel(wks, (/plot1,plot3/), (/1,2/), res_panel)
            else
                plot1 = gsn_csm_contour_map(wks, f->u(0,:,:), res1)
                plot2 = gsn_csm_contour_map(wks, f->v(0,:,:), res2)
                plot3 = gsn_csm_contour_map(wks, gh(0,:,:), res3)
                gsn_panel(wks, (/plot1,plot2,plot3/), (/2,2/), res_panel)
            end if
        end if
    end do

end
