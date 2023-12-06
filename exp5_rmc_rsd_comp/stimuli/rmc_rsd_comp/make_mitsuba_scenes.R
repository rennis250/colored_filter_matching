basen <- "mitsuba_comp_rmc_rsd"

illumns <- c("blue", "yellow", "red", "green")
bkgds <- c("voronoi_blue.png", "voronoi_yellow.png", "voronoi_red.png", "voronoi_green.png", "voronoi_diff_dist.png")

blue <- read.table("../../../base_stimuli/spectra/munsell_blue_lightest.spd", sep=' ')
yellow <- read.table("../../../base_stimuli/spectra/munsell_yellow_lightest.spd", sep=' ')
red <- read.table("../../../base_stimuli/spectra/munsell_red_lightest.spd", sep=' ')
green <- read.table("../../../base_stimuli/spectra/munsell_green_lightest.spd", sep=' ')

illumc <- 0
for(illumn in illumns) {
	illumc <- illumc + 1

	if(illumn == "blue") {
		illum <- red
    	illum[,2] <- blue[,2]
    	illumname <- paste0("./illums/blue.spd")
    	write.table(illum, illumname, quote=FALSE, sep=' ', row.names=FALSE, col.names=FALSE)
	} else if (illumn == "yellow") {
		illum <- red
    	illum[,2] <- yellow[,2]
    	illumname <- paste0("./illums/yellow.spd")
    	write.table(illum, illumname, quote=FALSE, sep=' ', row.names=FALSE, col.names=FALSE)
    } else if (illumn == "red") {
		illum <- red
    	illum[,2] <- red[,2]
    	illumname <- paste0("./illums/red.spd")
    	write.table(illum, illumname, quote=FALSE, sep=' ', row.names=FALSE, col.names=FALSE)
    } else if (illumn == "green") {
		illum <- red
    	illum[,2] <- green[,2]
    	illumname <- paste0("./illums/green.spd")
    	write.table(illum, illumname, quote=FALSE, sep=' ', row.names=FALSE, col.names=FALSE)
    }

	bkgdc <- 0
    for(bkgd in bkgds) {
    	bkgdc <- bkgdc + 1

    	ldc <- 0
		for(ld_scale in c(0.68, 1.0)) {
			ldc <- ldc + 1

			byc <- 0
		    for(by_mix in seq(-1,1,0.6666666)) {
		    	byc <- byc + 1

		    	rgc <- 0
		    	for(rg_mix in seq(-1,1,0.6666666)) {
		    		rgc <- rgc + 1

		    		refl <- red
		    		refl[,2] <- 0.5*(rg_mix*red[,2] + (1 - rg_mix)*green[,2]) + 0.5*(by_mix*blue[,2] + (1 - by_mix)*yellow[,2])
		    		refl[,2] <- ld_scale*refl[,2]
		    		refln <- paste0("./refls/rg_", rg_mix, "_by_", by_mix, "_ld_", ld_scale, ".spd")
		    		write.table(refl, refln, quote=FALSE, sep=' ', row.names=FALSE, col.names=FALSE)

		            template <- paste0('<?xml version=\'1.0\' encoding=\'utf-8\'?>
		            <!--

		                Automatically converted from COLLADA

		            -->

		            <scene version="0.6.0">
		                <integrator type="mlt"/>

		                <sensor type="perspective">
		                    <float name="farClip" value="100"/>
		                    <float name="fov" value="24.11"/>
		                    <string name="fovAxis" value="x"/>
		                    <float name="nearClip" value="0.1"/>

		                    <transform name="toWorld">
		                        <scale x="-1"/>
		                        <lookat target="0.5, 3, 1.5" origin="6, 12, 5" up="0, 0, 1.9"/>
		                    </transform>

		                    <sampler type="independent">
		                        <integer name="sampleCount" value="1000"/>
		                    </sampler>

		                    <film type="hdrfilm">
		                        <boolean name="highQualityEdges" value="true"/>
		                        <string name="fileFormat" value="openexr"/>
		                        <string name="componentFormat" value="float16"/>
		                        <string name="pixelFormat" value="rgb"/>
		                        <boolean name="banner" value="false"/>
		                        <integer name="height" value="400"/>
		                        <integer name="width" value="400"/>
		                        <rfilter type="gaussian">
		                            <float name="stddev" value="0.5"/>
		                        </rfilter>
		                    </film>
		                </sensor>

		                <bsdf type="diffuse" id="white">
		                    <spectrum name="reflectance" filename="../spectra/whiteWallspectrum.spd"/>
		                </bsdf>

		                <bsdf type="diffuse" id="colored_noise">
		                    <texture type="bitmap" name="reflectance">
		                        <string name="filename" value="./bkgds/', bkgd, '"/>
		                        <string name="wrapMode" value="mirror"/>
		                    </texture>
		                </bsdf>

		                <shape id="Plane-mesh_0" type="serialized">
		                    <ref id="colored_noise"/>
		                    <string name="filename" value="../blender/curveyBoxOverhead.serialized"/>
		                    <integer name="shapeIndex" value="0"/>
		                    <transform name="toWorld">
		                        <matrix value="5 0 0 0 0 5.5 0 0.5 0 0 5 0 0 0 0 1"/>
		                    </transform>
		                </shape>

		                <shape id="Plane_001-mesh_0" type="serialized">
		                    <ref id="colored_noise"/>
		                    <string name="filename" value="../blender/curveyBoxOverhead.serialized"/>
		                    <integer name="shapeIndex" value="1"/>
		                    <transform name="toWorld">
		                        <rotate x="1" angle="180"/>
		                        <matrix value="-2.18557e-07 0 5 5 0 5 0 0 -5 0 -2.18557e-07 5 0 0 0 1"/>
		                    </transform>
		                </shape>

		                <shape id="Plane_002-mesh_0" type="serialized">
		                    <ref id="colored_noise"/>
		                    <string name="filename" value="../blender/curveyBoxOverhead.serialized"/>
		                    <integer name="shapeIndex" value="2"/>
		                    <transform name="toWorld">
		                        <matrix value="-2.18557e-07 0 5 -5 0 5 0 0 -5 0 -2.18557e-07 5 0 0 0 1"/>
		                    </transform>
		                </shape>

		                <shape id="Plane_003-mesh_0" type="serialized">
		                    <ref id="colored_noise"/>
		                    <string name="filename" value="../blender/curveyBoxOverhead.serialized"/>
		                    <integer name="shapeIndex" value="3"/>
		                    <transform name="toWorld">
		                        <rotate x="1" angle="180"/>
		                        <matrix value="5 0 0 0 0 -2.18557e-07 -5 -5 0 5 -2.18557e-07 5 0 0 0 1"/>
		                    </transform>
		                </shape>

		    			<shape type="sphere">
		    				<transform name="toWorld">
		    					<scale value="0.2"/>
		    					<translate x="3.3" y="3.5" z="3.5"/>
		    				</transform>
		    				<emitter type="area">
		    					<spectrum name="radiance" filename="', illumname, '"/>
		    				</emitter>
		    			</shape>

		                <shape id="BigGlaven1_002-mesh_0" type="serialized">
		                    <bsdf type="dielectric">
		                        <spectrum name="specularTransmittance" filename="', refln, '"/>
		                        <spectrum name="specularReflectance" filename="', refln, '"/>
		                    </bsdf>
		                    <string name="filename" value="../blender/curveyBoxOverhead.serialized"/>
		                    <integer name="shapeIndex" value="6"/>
		                    <transform name="toWorld">
		    					<rotate x="1" angle="-90"/>
		    					<matrix value="0.100000000000000 0 0 0 0 0.100000000000000 0 3 0 0 0.100000000000000 1.5 0 0 0 1.000000000000000"/>
		                    </transform>
		                </shape>
		            </scene>')

		            fileConn <- file(paste0("scenes/", basen, "_illum_", as.character(illumc), "_rg_", as.character(rgc), "_by_", as.character(byc), "_ld_", as.character(ldc), "_bkgd_", substr(bkgd, 1, nchar(bkgd)-4), ".xml"))
		            writeLines(template, fileConn)
		            close(fileConn)
		        }
		    }
		}
	}
}
