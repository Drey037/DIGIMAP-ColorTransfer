var RBF_param_coff = 5;
function pinv(A)
{
    var z = numeric.svd(A), foo = z.S[0];
    var U = z.U, S = z.S, V = z.V;

    var m = A.length, n = A[0].length, tol = Math.max(m,n)*numeric.epsilon*foo,M = S.length;
    var i,Sinv = new Array(M);
    for (i=M-1;i!==-1;i--) { if(S[i]>tol) Sinv[i] = 1/S[i]; else Sinv[i] = 0; }
    return numeric.dot(numeric.dot(V,numeric.diag(Sinv)),numeric.transpose(U));
}
function solveEquation(A,b)
{
	var z = numeric.svd(A), foo = z.S[0];
	var U = z.U, S = z.S, V = z.V;

	var m = A.length, n = A[0].length, tol = Math.max(m,n)*numeric.epsilon*foo,M = S.length;
	var i,Sinv = new Array(M);
	for (i=M-1;i!==-1;i--) { if(S[i]>tol) Sinv[i] = 1/S[i]; else Sinv[i] = 0; }
	return numeric.dotMV( numeric.dot(numeric.dot(V,numeric.diag(Sinv)),numeric.transpose(U)), b );
}
function RGB2LAB(Q)
{
	var var_R = ( Q[0] / 255 );        //R from 0 to 255
	var var_G = ( Q[1] / 255 );        //G from 0 to 255
	var var_B = ( Q[2] / 255 );        //B from 0 to 255
	

	if ( var_R > 0.04045 ) var_R = Math.pow( ( var_R + 0.055 ) / 1.055, 2.4 );
	else                   var_R = var_R / 12.92;
	
	if ( var_G > 0.04045 ) var_G = Math.pow( ( var_G + 0.055 ) / 1.055, 2.4 );
	else                   var_G = var_G / 12.92;
	
	if ( var_B > 0.04045 ) var_B = Math.pow( ( var_B + 0.055 ) / 1.055, 2.4 );
	else                   var_B = var_B / 12.92;

	var_R = var_R * 100;
	var_G = var_G * 100;
	var_B = var_B * 100;

	var X = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805;
	var Y = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722;
	var Z = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505;
	
	var var_X = X / 95.047; 
	var var_Y = Y / 100; 
	var var_Z = Z / 108.883;

	if ( var_X > 0.008856 ) var_X = Math.pow(var_X, 1/3 );
	else                    var_X = ( 7.787 * var_X ) + ( 16 / 116 );
	
	if ( var_Y > 0.008856 ) var_Y = Math.pow(var_Y, 1/3 );
	else                    var_Y = ( 7.787 * var_Y ) + ( 16 / 116 );
	
	if ( var_Z > 0.008856 ) var_Z = Math.pow(var_Z, 1/3 );
	else                    var_Z = ( 7.787 * var_Z ) + ( 16 / 116 );

	var L = ( 116 * var_Y ) - 16;
	var A = 500 * ( var_X - var_Y );
	var B = 200 * ( var_Y - var_Z );
	
	return [L,A,B];
}

//LAB->XYZ->RGB from easyRGB.com 
function LAB2RGB(Q)
{
	//alert("L2R "+Q);
	var var_Y = ( Q[0] + 16 ) / 116;
	var var_X = Q[1] / 500 + var_Y;
	var var_Z = var_Y - Q[2] / 200;

	if ( var_Y > 0.206893034422 ) var_Y = Math.pow(var_Y,3);
	else                      var_Y = ( var_Y - 16 / 116 ) / 7.787;
	if ( var_X > 0.206893034422 ) var_X = Math.pow(var_X,3);
	else                      var_X = ( var_X - 16 / 116 ) / 7.787;
	if ( var_Z > 0.206893034422 ) var_Z = Math.pow(var_Z,3);
	else                      var_Z = ( var_Z - 16 / 116 ) / 7.787;

	var X = 95.047 * var_X;
	var Y = 100 * var_Y;
	var Z = 108.883 * var_Z;
	
	var_X = X / 100;
	var_Y = Y / 100;
	var_Z = Z / 100;

	var var_R = var_X *  3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
	var var_G = var_X * -0.9689 + var_Y *  1.8758 + var_Z *  0.0415;
	var var_B = var_X *  0.0557 + var_Y * -0.2040 + var_Z *  1.0570;

	if ( var_R > 0.0031308 ) var_R = 1.055 * Math.pow(var_R,1 / 2.4) - 0.055;
	else                     var_R = 12.92 * var_R;
	if ( var_G > 0.0031308 ) var_G = 1.055 * Math.pow(var_G,1 / 2.4) - 0.055;
	else                     var_G = 12.92 * var_G;
	if ( var_B > 0.0031308 ) var_B = 1.055 * Math.pow(var_B,1 / 2.4) - 0.055;
	else                     var_B = 12.92 * var_B;

	var R = var_R * 255;
	var G = var_G * 255;
	var B = var_B * 255;
	return [R,G,B];
}

/************************************************
* main function: rbfcolor
* param: 
*	image_name: image name
*	maxsize: 	max height/width of the processing image 
*	controlmap: used in the gallery version
*	callback: 	callback function
*	ngrid: 		grid size, n*n*n, default = 10 
*************************************************/
function rbfcolor(image_name, maxsize , controlmap, callback, grid_param){
	this.sum_grid_time = 0;
	this.sum_interp_time = 0;
	this.sum_draw_time = 0;
	
	var ngrid = grid_param || 12;
	
	// calculate the grid index and weight for and RGB color
	this.calculate_singlegrid = function(ori_r,ori_g,ori_b)
	{
        var ntmp = ngrid + 1;
        var ntmpsqr = ntmp*ntmp;
        var diff_x, diff_y,diff_z;
        var corner_ind;
        var tmpx, tmpy, tmpz;	
		
        tmpx = ori_r / 255 * ngrid;
        diff_x = tmpx - Math.floor(tmpx);
        tmpx = Math.floor(tmpx);
        if (tmpx == ngrid){
          tmpx = ngrid - 1;
          diff_x = 1;
        }
        tmpy = ori_g / 255 *ngrid;
        diff_y = tmpy - Math.floor(tmpy);
	  	tmpy = Math.floor(tmpy);
        if (tmpy == ngrid){
            tmpy = ngrid - 1;
            diff_y = 1;
        }
        tmpz = ori_b / 255 *ngrid;
        diff_z = tmpz - Math.floor(tmpz);
	  	tmpz = Math.floor(tmpz);
        if (tmpz == ngrid){
            tmpz = ngrid - 1;
            diff_z = 1;
        }
   
        corner_ind = tmpz * ntmpsqr + tmpy * ntmp +tmpx;

		var res = new Array(16);
		
		res[0] = corner_ind;
        res[1] = corner_ind+ntmpsqr;
        res[2] = corner_ind+ntmp;
        res[3] = corner_ind+ntmp+ntmpsqr;
        res[4] = corner_ind+1;
        res[5] = corner_ind+ntmpsqr+1;
        res[6] = corner_ind+ntmp+1;
        res[7] = corner_ind+ntmp+ntmpsqr+1;
 
        res[8] = (1-diff_x)*(1-diff_y)*(1-diff_z);
        res[9] = (1-diff_x)*(1-diff_y)*diff_z;
        res[10] = (1-diff_x)*diff_y*(1-diff_z);
        res[11] = (1-diff_x)*diff_y*diff_z;
        res[12] = diff_x*(1-diff_y)*(1-diff_z);
        res[13] = diff_x*(1-diff_y)*diff_z;
        res[14] = diff_x*diff_y*(1-diff_z);
        res[15] = diff_x*diff_y*diff_z;
		
		return res;
	}
	// preprocessing grid for an input
	this.prepare_grid = function()
	{
		// initialize the LAB grid
		this.grid_size = (ngrid+1) * (ngrid+1) * (ngrid+1);
		this.grid_lab = new Array(this.grid_size);
		
		var step = 255.0 / ngrid;
		var tot = 0;
		for(var i=0;i<ngrid+1;i++)
		for(var j=0;j<ngrid+1;j++)
		for(var k=0;k<ngrid+1;k++)
		{
			this.grid_lab[tot] = RGB2LAB([k*step,j*step,i*step]);
		    tot = tot + 1;
		}
		
		// calculate the weight index and map for image pixel
        this.weightindex = new Array(this.img_area);
        this.weightmap   = new Array(this.img_area);

        for (var i=0; i<this.img_area; i++) {
			var res = this.calculate_singlegrid( this.ori_r[i],this.ori_g[i],this.ori_b[i]);
      
  		  this.weightindex[i] = new Array(8);
		  for(var j=0;j < 8;j++)
		  {
		  	this.weightindex[i][j] = res[j];
		  }
		  this.weightmap[i] = new Array(8);
		  for(var j=0;j < 8;j++)
		  {
		  	this.weightmap[i][j] = res[j+8];
		  }
        }
		// initialize result grid
		this.grid_R = new Array(this.grid_size);
		this.grid_G = new Array(this.grid_size);
		this.grid_B = new Array(this.grid_size);
	}
	
	var imageObj = new Image();
	var ptr = this;
	
	imageObj.onload = function(){
		ptr.originalWidth = imageObj.width;
		ptr.originalHeight = imageObj.height;
		
		// paint on a virtual canvas with maxsize 
		var img_width, img_height;
		img_height = imageObj.height;
		img_width = imageObj.width;
		
		

		ptr.img_width = img_width;
		ptr.img_height = img_height;

		var newCanvas = $("<canvas>")[0];
		newCanvas.width = img_width;
		newCanvas.height = img_height;
		
		var newContext = newCanvas.getContext("2d");
		newContext.width = img_width;
		newContext.height = img_height;
		
        newContext.drawImage(imageObj, 0,0,img_width,img_height);
		
        var imgData = newContext.getImageData(0,0,img_width,img_height);
        ptr.img_area = img_width * img_height;

        ptr.ori_b = new Array(ptr.img_area);
        ptr.ori_g = new Array(ptr.img_area);
        ptr.ori_r = new Array(ptr.img_area);
        var cnt = 0;
        for (var i=0; i<imgData.data.length; i+=4)
        {
            ptr.ori_r[cnt] = imgData.data[i];
            ptr.ori_g[cnt] = imgData.data[i+1];
            ptr.ori_b[cnt] = imgData.data[i+2];
            cnt = cnt + 1;
        }
		
		ptr.prepare_grid();
		
        ptr.res_b = new Array(ptr.img_area);
        ptr.res_g = new Array(ptr.img_area);
        ptr.res_r = new Array(ptr.img_area);
		
		
		ptr.mask = new Array(ptr.img_area);
		for(var i=0;i < ptr.img_area;i++)
		{
			ptr.mask[i] = 1.0;
		}
		callback();
	};
	imageObj.src = image_name;
	
	
	this.updateMask = function(image_name)
	{
		var imageObj = new Image();
		var ptr = this;
	
		imageObj.onload = function(){
			// paint on a virtual canvas
			var img_width, img_height;

			img_width = ptr.img_width;
			img_height = ptr.img_height;

			var newCanvas = $("<canvas>")[0];
			newCanvas.width = img_width;
			newCanvas.height = img_height;
		
			var newContext = newCanvas.getContext("2d");
			newContext.width = img_width;
			newContext.height = img_height;
		
	        newContext.drawImage(imageObj, 0,0,img_width,img_height);
		
	        var imgData = newContext.getImageData(0,0,img_width,img_height);

	        var cnt = 0;
	        for (var i=0; i<imgData.data.length; i+=4)
	        {
	            ptr.mask[cnt] = Math.max( Math.max( imgData.data[i],imgData.data[i+1] ), imgData.data[i+2] ) / 255;
	            cnt = cnt + 1;
	        }
		
			alert("mask loaded");
		};
		imageObj.src = image_name;
	}
	

	this.CalculateLABDistance = function(l1,a1,b1,l2,a2,b2)
	{
		
		var K1 = 0.045, K2 = 0.015;
		var del_L = l1 - l2;
		var c1 = Math.sqrt(a1*a1+b1*b1); 
		var c2 = Math.sqrt(a2*a2+b2*b2);
		var c_ab = c1 - c2;
		var h_ab = (a1-a2)*(a1-a2)+(b1-b2)*(b1-b2) - c_ab*c_ab;
		return del_L*del_L + c_ab *c_ab /(1+K1*c1)/(1+K1*c1) + h_ab / (1+K2*c1)/(1+K2*c1);
		
	}
	// variant RBF interpolation
	this.CalculateSinglePoint = function(param, matrixsize, oldpalette_L,oldpalette_A,oldpalette_B, diffpalette_L, diffpalette_A, diffpalette_B, vsrc )
	{
		var tmpMat = new Array(matrixsize);
		for(var i=0;i < matrixsize;i++ )
		{
			tmpMat[i] = new Array(matrixsize);
		}
		for(var u=0;u < matrixsize;u++)
		for(var v=0;v < matrixsize;v++)
		{
			var r = this.CalculateLABDistance(oldpalette_L[u],oldpalette_A[u],oldpalette_B[u],oldpalette_L[v],oldpalette_A[v],oldpalette_B[v]);		   
		    tmpMat[u][v] = Math.exp(-r*param);
		}
		var tmpD = new Array(matrixsize);
		for(var u=0;u < matrixsize;u++)
		{
			var r = this.CalculateLABDistance(oldpalette_L[u],oldpalette_A[u],oldpalette_B[u],vsrc[0],vsrc[1],vsrc[2]);
			tmpD[u] = Math.exp(-r*param);
		}

		var precompute_pinv = numeric.dotMV( pinv(tmpMat), tmpD );

		var delta_L = 0;
		var delta_A = 0;
		var delta_B = 0;

		var scale = 0;
		for(var j = 0;j < matrixsize;j++ )
		{
			scale = scale + Math.max( precompute_pinv[j], 0.0 );
		}
		for( var j = 0;j < matrixsize;j++ ) if(precompute_pinv[j]>0)
		{
			delta_L = delta_L + precompute_pinv[j] / scale * diffpalette_L[j];
			delta_A = delta_A + precompute_pinv[j] / scale * diffpalette_A[j];
			delta_B = delta_B + precompute_pinv[j] / scale * diffpalette_B[j];
		}
		
		return [ vsrc[0] + delta_L , vsrc[1] + delta_A, vsrc[2]+delta_B ];
	}

	this.Bigchange = function(dir)
	{
		return ( Math.abs(dir[0])>0.5 || Math.abs(dir[1])>0.5 || Math.abs(dir[2])>0.5 );
	}
	
	this.OutBoundary = function(testrgb )
	{
		var out_threshold = 0.5;
		return ( testrgb[0] < -out_threshold || testrgb[0] > 255+out_threshold ||
		    testrgb[1] < -out_threshold || testrgb[1] > 255+out_threshold ||
		    testrgb[2] < -out_threshold || testrgb[2] > 255+out_threshold  );
	}
	
	this.FindBoundary = function(vsrc, dir, l, r)
	{
		//Assume dir is large
		var mid;
		for(var iter = 0;iter < 15;iter++)
		{
			mid = 0.5 * (l+r);
			var testrgb = LAB2RGB([ vsrc[0] + mid * dir[0], vsrc[1] + mid * dir[1], vsrc[2] + mid * dir[2] ]);
			
			if(this.OutBoundary(testrgb))
			{
				r = mid;
			}else
			{
				l = mid;
			}
		}
		return l;
	}

	this.CalculateGridResult = function(palette_size, oldpalette_L,oldpalette_A,oldpalette_B, diffpalette_L, diffpalette_A, diffpalette_B )
	{
		//average distance
		var tot = 0;
		var totr = 0;
		for(var u=0;u < palette_size;u++)
		for(var v=u+1;v < palette_size;v++)
		{
			var r = this.CalculateLABDistance(oldpalette_L[u],oldpalette_A[u],oldpalette_B[u],oldpalette_L[v],oldpalette_A[v],oldpalette_B[v]);
			
			tot = tot + 1;
			totr = totr + Math.sqrt(r);
		}
		
		if( palette_size > 1 )
		{
			totr = totr / tot;
		}else
		{
			totr = 1.0;
		}
		
		var param = RBF_param_coff / (totr*totr);
		
		for(var i=0;i < this.grid_size;i++ )
		{
			var vsrc = this.grid_lab[i];
			
			var tdiff_L = new Array(palette_size);
			var tdiff_A = new Array(palette_size);
			var tdiff_B = new Array(palette_size);
			
			for(var j=0; j < palette_size;j++)
			{
				
				var dir = [ diffpalette_L[j],diffpalette_A[j],diffpalette_B[j] ];

				if(this.Bigchange(dir))
				{
					var pc = this.LAB2RGB( [vsrc[0]+dir[0],vsrc[1]+dir[1],vsrc[2]+dir[2]] );
					if(this.OutBoundary(pc))
					{
						var M = [ oldpalette_L[j]+dir[0], oldpalette_A[j]+dir[1], oldpalette_B[j]+dir[2] ];
						var Mdir = [vsrc[0]-oldpalette_L[j],vsrc[1]-oldpalette_A[j],vsrc[2]-oldpalette_B[j]];
						
						var t1 = this.FindBoundary( M , Mdir,0,1 );
						
						if(!this.OutBoundary(this.LAB2RGB(  [M[0]+Mdir[0],M[1]+Mdir[1],M[2]+Mdir[2]]  )))
						{
							alert("M+Mdir Error");
						}
						
						
						var t2 = this.FindBoundary( [ oldpalette_L[j], oldpalette_A[j], oldpalette_B[j] ], dir, 1,300 );	
						
						tdiff_L[j] = (dir[0] - (1-t1)*Mdir[0]) ;/// t2;
						tdiff_A[j] = (dir[1] - (1-t1)*Mdir[1]) / t2;
						tdiff_B[j] = (dir[2] - (1-t1)*Mdir[2]) / t2;
						
						if(t1>1)
						{
							alert("t1>1 at case 1 "+t1+" "+Mdir);	
						}
						if(t2<1)
						{
							alert("t2<1 at case 1 "+t2);
						}
					}else
					{
						var t1 = this.FindBoundary( vsrc, dir,1,300 );
						var t2 = this.FindBoundary( [ oldpalette_L[j], oldpalette_A[j], oldpalette_B[j] ], dir, 1,300 );	
						
						if(t2<1)
						{
							alert("t2<1 at case 2 "+t2);
						}
						var lambda = Math.min( t1/t2 , 1.0 );
						tdiff_L[j] = diffpalette_L[j] ;//* lambda;
						tdiff_A[j] = diffpalette_A[j] * lambda;
						tdiff_B[j] = diffpalette_B[j] * lambda;
					}
				}else
				{
					tdiff_L[j] = diffpalette_L[j];
					tdiff_A[j] = diffpalette_A[j];
					tdiff_B[j] = diffpalette_B[j];
				}
				
			}
			
			var res = this.CalculateSinglePoint( param, palette_size, oldpalette_L,oldpalette_A,oldpalette_B, tdiff_L, tdiff_A, tdiff_B, vsrc );
			
			var pc = LAB2RGB(res);
			this.grid_R[i] = pc[0];
			this.grid_G[i] = pc[1];
			this.grid_B[i] = pc[2];			
		
			this.grid_R[i] = Math.max(0, Math.min(this.grid_R[i],255))
			this.grid_G[i] = Math.max(0, Math.min(this.grid_G[i],255))
			this.grid_B[i] = Math.max(0, Math.min(this.grid_B[i],255))
		}
	}
	
	this.drawResImage = function(output_canvas)
	{
		var output_context = output_canvas.getContext('2d');
        var output_imgData = output_context.createImageData(this.img_width, this.img_height);
        var cnt = 0;
        for (var i=0; i<output_imgData.data.length; i+=4)
        {
            output_imgData.data[i] = this.res_r[cnt] * this.mask[cnt] + this.ori_r[cnt] * (1-this.mask[cnt]);
            output_imgData.data[i+1]= this.res_g[cnt] * this.mask[cnt] + this.ori_g[cnt] * (1-this.mask[cnt]);
            output_imgData.data[i+2]= this.res_b[cnt] * this.mask[cnt] + this.ori_b[cnt] * (1-this.mask[cnt]);
            output_imgData.data[i+3] = 255;
            cnt = cnt + 1;
        }
		var newCanvas = $("<canvas>")
        .attr("width", output_imgData.width)
        .attr("height", output_imgData.height)[0];
		
		newCanvas.getContext("2d").putImageData(output_imgData, 0, 0);
		output_context.scale(output_canvas.width / this.img_width,output_canvas.height / this.img_height);
		output_context.drawImage(newCanvas, 0, 0);
		output_context.scale(this.img_width / output_canvas.width , this.img_height / output_canvas.height );
		
	}
	this.drawRes = function(input_palette, output_palette, palette_size, output_canvas)
	{
		var oldpalette_L = new Array(palette_size+1);
		var oldpalette_A = new Array(palette_size+1);
		var oldpalette_B = new Array(palette_size+1);
		var newpalette_L = new Array(palette_size+1);
		var newpalette_A = new Array(palette_size+1);
		var newpalette_B = new Array(palette_size+1);
		
		var tmplab;
		for(var i=0;i<palette_size;i++)
		{
			tmplab = RGB2LAB( [ input_palette[i*3]*255,input_palette[i*3+1]*255,input_palette[i*3+2]*255 ] );
			oldpalette_L[i] = tmplab[0];oldpalette_A[i] = tmplab[1];oldpalette_B[i] = tmplab[2];
			
			tmplab = RGB2LAB( [ output_palette[i*3]*255,output_palette[i*3+1]*255,output_palette[i*3+2]*255 ] );
			newpalette_L[i] = tmplab[0];newpalette_A[i] = tmplab[1];newpalette_B[i] = tmplab[2];
		}
		
		oldpalette_L[palette_size]=0;
		oldpalette_A[palette_size]=0;
		oldpalette_B[palette_size]=0;
		
		var diffpalette_L = new Array(palette_size+1);
		var diffpalette_A = new Array(palette_size+1);
		var diffpalette_B = new Array(palette_size+1);
		
		for(var i=0;i < palette_size;i++)
		{
			diffpalette_L[i] = newpalette_L[i]-oldpalette_L[i];
			diffpalette_A[i] = newpalette_A[i]-oldpalette_A[i];
			diffpalette_B[i] = newpalette_B[i]-oldpalette_B[i];
		}
		
		diffpalette_L[palette_size]=0;
		diffpalette_A[palette_size]=0;
		diffpalette_B[palette_size]=0;
		
		this.CalculateGridResult(palette_size, oldpalette_L,oldpalette_A,oldpalette_B, diffpalette_L, diffpalette_A, diffpalette_B );
		
        for (var i=0; i<this.img_area;i++){
			
			var tmpR = 0;
			var tmpG = 0;
			var tmpB = 0;
         	for (var k=0; k<8; k++){
             	tmpR = tmpR + this.grid_R[this.weightindex[i][k]] * this.weightmap[i][k];
            	tmpG = tmpG + this.grid_G[this.weightindex[i][k]] * this.weightmap[i][k];
            	tmpB = tmpB + this.grid_B[this.weightindex[i][k]] * this.weightmap[i][k];
        	}
			this.res_r[i] = tmpR;
			this.res_g[i] = tmpG;
			this.res_b[i] = tmpB;
 		}
		
		this.drawResImage(output_canvas);
	}
	
	this.drawImage = function(output_canvas)
	{
		var output_context = output_canvas.getContext('2d');
        var output_imgData = output_context.createImageData(this.img_width, this.img_height);
        var cnt = 0;
        for (var i=0; i<output_imgData.data.length; i+=4)
        {
            output_imgData.data[i] = this.ori_r[cnt];
            output_imgData.data[i+1]= this.ori_g[cnt];
            output_imgData.data[i+2]= this.ori_b[cnt];
            output_imgData.data[i+3] = 255;
            cnt = cnt + 1;
        }
		var newCanvas = $("<canvas>")
        .attr("width", output_imgData.width)
        .attr("height", output_imgData.height)[0];
		
		newCanvas.getContext("2d").putImageData(output_imgData, 0, 0);
		
		output_context.scale(output_canvas.width / this.img_width,output_canvas.height / this.img_height);
		output_context.drawImage(newCanvas, 0, 0);
		output_context.scale(this.img_width / output_canvas.width , this.img_height / output_canvas.height );
	}
	
	
	//other function
	this.RGB2LAB = function(Q)
	{
		return RGB2LAB(Q);
	}
	this.LAB2RGB = function(Q)
	{
		return LAB2RGB(Q);
	}

	
	
	//Kmeans
	this.weight_gridacc_kmeans = function(center_num, weights) {
		var ngrid = 16;
		var grid_size = ngrid * ngrid * ngrid;
		var step_size = 255.0 / (ngrid - 1);
		var sample_cnt = new Array(grid_size);
		var sample_sum = new Array(grid_size);
		var sample_weight_sum = new Array(grid_size); // New array to store the sum of weights
	
		for (var i = 0; i < grid_size; i++) {
			sample_cnt[i] = 0;
			sample_sum[i] = new Array(3);
			sample_sum[i][0] = 0;
			sample_sum[i][1] = 0;
			sample_sum[i][2] = 0;
			sample_weight_sum[i] = 0; // Initialize weight sum
		}
	
		for (var i = 0; i < this.img_area; i++) {
			var p = RGB2LAB([this.ori_r[i], this.ori_g[i], this.ori_b[i]]);
			var bin1 = Math.round(this.ori_r[i] / step_size);
			var bin2 = Math.round(this.ori_g[i] / step_size);
			var bin3 = Math.round(this.ori_b[i] / step_size);
			var bin = bin1 * ngrid * ngrid + bin2 * ngrid + bin3;
	
			sample_cnt[bin] += 1;
			sample_sum[bin][0] += p[0];
			sample_sum[bin][1] += p[1];
			sample_sum[bin][2] += p[2];
			sample_weight_sum[bin] += weights[i]; // Add the weight of the sample
		}
	
		var tot = 0;
		for (var i = 0; i < grid_size; i++) {
			if (sample_cnt[i] > 0) {
				tot++;
			}
		}
	
		var D = new Array(tot);
		tot = 0;
		for (var i = 0; i < grid_size; i++) {
			if (sample_cnt[i] > 0) {
				D[tot] = new Array(3);
				D[tot][0] = sample_sum[i][0] / sample_cnt[i];
				D[tot][1] = sample_sum[i][1] / sample_cnt[i];
				D[tot][2] = sample_sum[i][2] / sample_cnt[i];
				sample_cnt[tot] = sample_cnt[i];
				tot++;
			}
		}
	
		var Center = new Array(center_num + 1);
		var pickcnt = new Array(tot);
		for (var i = 0; i < tot; i++) {
			pickcnt[i] = sample_cnt[i];
		}
	
		for (var i = 0; i < center_num; i++) {
			var idx = 0;
			for (var j = 0; j < tot; j++) {
				if (pickcnt[j] > pickcnt[idx]) {
					idx = j;
				}
			}
	
			Center[i] = D[idx];
			for (var j = 0; j < tot; j++) {
				var dis = 0;
				for (var k = 0; k < 3; k++) {
					dis += (D[idx][k] - D[j][k]) * (D[idx][k] - D[j][k]);
				}
				dis /= (80 * 80);
	
				// Adjust the weight update based on the sample's weight
				pickcnt[j] *= (1 - Math.exp(-dis)) * (sample_weight_sum[j] / sample_cnt[j]);
			}
		}
	
		// Add black
		Center[center_num] = RGB2LAB([0, 0, 0]);
		center_num = center_num + 1;
	
		var cnt = new Array(center_num);
		var sumD = new Array(center_num);
	
		for (var iter = 0; iter < 20; iter++) {
			for (var i = 0; i < center_num; i++) {
				sumD[i] = [0, 0, 0];
				cnt[i] = 0;
			}
	
			for (var i = 0; i < tot; i++) {
				var min_id = -1;
				var min_v = 1e100;
				for (var j = 0; j < center_num; j++) {
					var r = (D[i][0] - Center[j][0]) * (D[i][0] - Center[j][0]) + (D[i][1] - Center[j][1]) * (D[i][1] - Center[j][1]) + (D[i][2] - Center[j][2]) * (D[i][2] - Center[j][2]);
					if (r < min_v) {
						min_v = r;
						min_id = j;
					}
				}
	
				cnt[min_id] += sample_cnt[i];
				sumD[min_id][0] += sample_cnt[i] * D[i][0];
				sumD[min_id][1] += sample_cnt[i] * D[i][1];
				sumD[min_id][2] += sample_cnt[i] * D[i][2];
			}
	
			for (var i = 0; i < center_num; i++) {
				if (cnt[i] > 0) {
					Center[i][0] = sumD[i][0] / cnt[i];
					Center[i][1] = sumD[i][1] / cnt[i];
					Center[i][2] = sumD[i][2] / cnt[i];
				}
			}
		}
	
		center_num = center_num - 1;
		var res = new Array(3 * center_num);
		for (var i = 0; i < center_num; i++) {
			var ps = LAB2RGB(Center[i]);
			res[i * 3 + 0] = Math.max(0.0, Math.min(1.0, ps[0] / 255));
			res[i * 3 + 1] = Math.max(0.0, Math.min(1.0, ps[1] / 255));
			res[i * 3 + 2] = Math.max(0.0, Math.min(1.0, ps[2] / 255));
		}
	
		return res;
	}

	this.gridacc_kmeans = function(center_num)
	{
          var ngrid = 16;

          var grid_size = ngrid * ngrid * ngrid;

          var step_size = 255.0 / (ngrid-1);
          var sample_cnt = new Array(grid_size);
          var sample_sum = new Array(grid_size);
          for(var i=0;i < grid_size;i++)
          {
              sample_cnt[i] = 0;
              sample_sum[i] = new Array(3);
              sample_sum[i][0] = 0;
              sample_sum[i][1] = 0;
              sample_sum[i][2] = 0;
          }

          for( var i = 0;i < this.img_area;i++ )
          {
              var p = RGB2LAB([ this.ori_r[i], this.ori_g[i], this.ori_b[i] ]);

              var bin1 = Math.round( this.ori_r[i] / step_size );
              var bin2 = Math.round( this.ori_g[i] / step_size );
              var bin3 = Math.round( this.ori_b[i] / step_size );

              var bin = bin1 * ngrid * ngrid + bin2 * ngrid + bin3;

              sample_cnt[bin] = sample_cnt[bin] + 1;
              sample_sum[bin][0] = sample_sum[bin][0] + p[0];
              sample_sum[bin][1] = sample_sum[bin][1] + p[1];
              sample_sum[bin][2] = sample_sum[bin][2] + p[2];
          }

          var tot = 0;
          for(var i=0;i < grid_size;i++) if(sample_cnt[i]>0)
          {
              tot = tot + 1;
          }

          var D = new Array(tot);
          tot = 0;
          for(var i=0;i < grid_size;i++) if(sample_cnt[i]>0)
          {
              D[tot] = new Array(3);
              D[tot][0] = sample_sum[i][0] / sample_cnt[i];
              D[tot][1] = sample_sum[i][1] / sample_cnt[i];
              D[tot][2] = sample_sum[i][2] / sample_cnt[i];
              sample_cnt[tot] = sample_cnt[i];

              tot = tot + 1;
          }

          var Center = new Array(center_num+1);

          var pickcnt = new Array(tot);
          for(var i=0;i < tot;i++)
          {
              pickcnt[i] = sample_cnt[i];
          }

          for (var i = 0;i < center_num;i++)
          {
              var idx = 0;
              for(var j=0;j<tot;j++) if(pickcnt[j] > pickcnt[idx])
              {
                  idx = j;
              }

              Center[i] = D[idx];
              for(var j=0; j < tot;j++)
              {
                  var dis = 0;
                  for(var k=0;k<3;k++)
                  {
                      dis = dis + (D[idx][k] - D[j][k])*(D[idx][k] - D[j][k]);
                  }
                  dis = dis / (80 * 80);

                  pickcnt[j] = pickcnt[j] * ( 1 - Math.exp(-dis) );
              }
          }

          //add black
          Center[center_num] = RGB2LAB([0,0,0]);
          center_num = center_num + 1;


          var cnt = new Array(center_num);
          var sumD = new Array(center_num);

          for( var iter=0;iter < 20;iter++)
          {
              for(var i=0;i < center_num;i++)
              {
                  sumD[i] = [0,0,0];
                  cnt[i] = 0;
              }

              for(var i = 0;i < tot;i++)
              {
                  var min_id = -1;
                  var min_v = 1e100;
                  for( var j = 0; j < center_num;j++ )
                  {
                      var r = (D[i][0]-Center[j][0])*(D[i][0]-Center[j][0])+(D[i][1]-Center[j][1])*(D[i][1]-Center[j][1])+(D[i][2]-Center[j][2])*(D[i][2]-Center[j][2]);
                      if( r < min_v )
                      {
                          min_v = r;
                          min_id = j;
                      }
                  }

                  cnt[min_id] = cnt[min_id] + sample_cnt[i];
                  sumD[min_id][0] = sumD[min_id][0] + sample_cnt[i] * D[i][0];
                  sumD[min_id][1] = sumD[min_id][1] + sample_cnt[i] * D[i][1];
                  sumD[min_id][2] = sumD[min_id][2] + sample_cnt[i] * D[i][2];
              }


              for( var i = 0;i < center_num;i++ ) if(cnt[i]>0)
              {
                  Center[i][0] = sumD[i][0] / cnt[i];
                  Center[i][1] = sumD[i][1] / cnt[i];
                  Center[i][2] = sumD[i][2] / cnt[i];
              }
          }


          center_num = center_num - 1;
          var res = new Array(3*center_num);
          for(var i = 0;i < center_num;i++ )
          {
              var ps = LAB2RGB(Center[i]);
              res[i*3+0] = Math.max(0.0, Math.min(1.0, ps[0]/255) );
              res[i*3+1] = Math.max(0.0, Math.min(1.0, ps[1]/255) );
              res[i*3+2] = Math.max(0.0, Math.min(1.0, ps[2]/255) );
          }

          return res;
      }
}