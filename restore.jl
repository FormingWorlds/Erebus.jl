                coords=coords,
                dsubgrids=dsubgrids,
            )

            # ---------------------------------------------------------------------
            # interpolate DSXX, DSXY to markers
            # ---------------------------------------------------------------------
            update_marker_stress!(xm, ym, sxxm, sxym, DSXX, DSXY, marknum; coords=coords)

            # ---------------------------------------------------------------------
            # apply subgrid temperature diffusion on markers,
            # compute DTsubgrid
            # ---------------------------------------------------------------------
            apply_subgrid_temperature_diffusion!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                tk1,
                DT,
                TKSUM,
                RHOCPSUM,
                dt,
                marknum,
                marker_property_mode;
                coords=coords,
                dsubgridt=dsubgridt,
            )

            # ---------------------------------------------------------------------
