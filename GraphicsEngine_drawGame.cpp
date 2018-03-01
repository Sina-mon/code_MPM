#include "GraphicsEngine.h"

// ----------------------------------------------------------------------------
void GraphicsEngine::drawGame(void)
{
	// global colors
	glm::vec4 f4Color_Default	= _BLUE;
	glm::vec4 f4Color_State_Mid	= _RED;
	glm::vec4 f4Color_State_Max	= _RED;
	glm::vec4 f4Color_DisplacementControl	= _GRAY;
	glm::vec4 f4Color_Surface	= _GREEN;
	glm::vec4 f4Color_ESO_False	= _GRAY;

	if(true)
	{// create shadow map
		gl_Shadow_Texture->bindRenderTarget();
//		gl_Canvas_Texture->bindRenderTarget(i_ScreenWidth, i_ScreenHeight);
		gl_ShadowProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		GLuint transformationShadowLocation = gl_ShadowProgram.getUniformLocation("transformationShadowMatrix");

		std::vector<MaterialPoint_CPDI_CC *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints_CPDI();
//		std::vector<MaterialPoint_BC *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();

		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];
			// particle position
			float fSize = glm::pow(thisMP->d_Volume, 1.0/3.0);
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection() * glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// draw
			gl_Particle_Mesh->Draw();
		}

		gl_ShadowProgram.unuse();
	}

	if(true)
	{// draw all objects to the canvas (not the screen)
		v_Canvas_Texture[(int)enum_Canvas::MAIN]->bindRenderTarget();
//		gl_Canvas_Texture->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Diffuse_Texture->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		// material points, classic -------------------------------------------
		std::vector<MaterialPoint_BC *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 1.0*(mpm_PhysicsEngine->d3_Length_World.x) / 1.0e3;
//			float fSize = 0.5*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			glm::vec4 f4objectColor = _RED;
			if(thisMP->b_Mark_Stress)
				f4objectColor = _GREEN;
			if(thisMP->b_Surface)
			{
				f4objectColor = _BLUE;
				fSize *= 1.0;
			}
			if(thisMP->b3_DisplacementControl.x || thisMP->b3_DisplacementControl.y || thisMP->b3_DisplacementControl.z)
				f4objectColor = _GRAY;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}
		// material points, CPDI ----------------------------------------------
		std::vector<MaterialPoint_CPDI_CC *> vMaterialPoint_CPDI = mpm_PhysicsEngine->getMaterialPoints_CPDI();
		for(int index_MP = 0; index_MP < vMaterialPoint_CPDI.size(); index_MP++)
		{
			MaterialPoint_CPDI_CC *thisMP = vMaterialPoint_CPDI[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			glm::vec4 f4objectColor = _BLUE;
			if(thisMP->b_Mark_Stress)
				f4objectColor = _RED;
			if(thisMP->b_Surface)
				f4objectColor = _BLUE;
			if(thisMP->b3_DisplacementControl.x || thisMP->b3_DisplacementControl.y || thisMP->b3_DisplacementControl.z)
				f4objectColor = _BLACK;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();

			for(int index_Corner = 0; index_Corner < 4; index_Corner++)
			{
				// particle position
				float fSize = 1.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
				// particle color
				//glm::vec4 f4objectColor = 0.5f*_BLUE + 0.5f*_WHITE;
				glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

				// shadow
				glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
				glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				// camera and model transformation matices
				Transformation glTransformation(thisMP->a_Corner[index_Corner].d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
				glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
				glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
				glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

				gl_Particle_Mesh->Draw();
			}
		}

		// grid points --------------------------------------------------------
		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();
		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(thisGP->b3_Fixed == glm::bvec3{false, false, false})
				continue;

			// particle position
			float fSize = 5.0*(mpm_PhysicsEngine->d3_Length_World.x) / 1.0e4;
			glm::vec3 f3Size = glm::vec3(fSize,fSize,fSize);
			Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), f3Size);
			// particle color
			glm::vec4 f4objectColor = _GRAY;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		gl_BasicProgram.unuse();
	}

	if((int)enum_Canvas::SOLID < (int)enum_Canvas::COUNT)
	{// Solid, only containing the material points
		v_Canvas_Texture[(int)enum_Canvas::SOLID]->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Diffuse_Texture->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<MaterialPoint_BC *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();
		std::vector<MaterialPoint_CPDI_CC *> vMaterialPoint_CPDI = mpm_PhysicsEngine->getMaterialPoints_CPDI();

		// material points, classic -------------------------------------------
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			glm::vec4 f4objectColor = _BLUE;
			if(thisMP->b_Mark_Stress)
				f4objectColor = _BLUE;
			if(thisMP->b3_DisplacementControl.x || thisMP->b3_DisplacementControl.y || thisMP->b3_DisplacementControl.z)
				f4objectColor = _GRAY;
			if(thisMP->b_Surface)
				f4objectColor = _GREEN;

			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}
		// material points, CPDI ----------------------------------------------
		for(int index_MP = 0; index_MP < vMaterialPoint_CPDI.size(); index_MP++)
		{
			MaterialPoint_CPDI_CC *thisMP = vMaterialPoint_CPDI[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			glm::vec4 f4objectColor = _BLUE;
			if(thisMP->b_Mark_Stress)
				f4objectColor = _BLUE;
			if(thisMP->b_Surface)
				f4objectColor = _GREEN;
			if(thisMP->b3_DisplacementControl.x || thisMP->b3_DisplacementControl.y || thisMP->b3_DisplacementControl.z)
				f4objectColor = _GRAY;

			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();

			for(int index_Corner = 0; index_Corner < 4; index_Corner++)
			{
				// particle position
				float fSize = 1.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);

				// shadow
				glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
				glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				// camera and model transformation matices
				Transformation glTransformation(thisMP->a_Corner[index_Corner].d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
				glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
				glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
				glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

				gl_Particle_Mesh->Draw();
			}
		}

		// grid points --------------------------------------------------------
		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();
		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(thisGP->b3_Fixed.y == false && thisGP->b3_Fixed.x == false)
				continue;

			// particle position
			float fSize = 5.0*(mpm_PhysicsEngine->d3_Length_World.x) / 1.0e4;
			glm::vec3 f3Size = glm::vec3(fSize,fSize,fSize);
			if(thisGP->b3_Fixed.y == true)
			{
				f3Size.x *= 10.0;
				f3Size.z *= 10.0;
			}
			else if(thisGP->b3_Fixed.x == true)
			{
				f3Size.y *= 10.0;
				f3Size.z *= 10.0;
			}
			Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), f3Size);
			// particle color
			glm::vec4 f4objectColor = _BLACK;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		gl_BasicProgram.unuse();
	}

	if((int)enum_Canvas::J2_PLASTICSTRAIN < (int)enum_Canvas::COUNT)
	{// MP parameteric value, equivalent plastic strain
		v_Canvas_Texture[(int)enum_Canvas::J2_PLASTICSTRAIN]->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Diffuse_Texture->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<MaterialPoint_BC *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();
		std::vector<MaterialPoint_CPDI_CC *> vMaterialPoint_CPDI = mpm_PhysicsEngine->getMaterialPoints_CPDI();

//		vMaterialPoint.insert(vMaterialPoint.end(), vMaterialPoint_CPDI.begin(), vMaterialPoint_CPDI.end());
		// material points ----------------------------------------------------
		ConstitutiveRelation CR;
//		double d6Strain_Minimum[6] = {0, 0, 0, 0, 0, 0};
//		double d6Strain_Maximum[6] = {0, 0, 0, 0, 0, 0};
//		if(vMaterialPoint.size() != 0)
//		{
//			d6Strain_Minimum[0] = 0.0;//vMaterialPoint[0]->d_YieldStress/vMaterialPoint[0]->d_ElasticModulus;
//			d6Strain_Maximum[0] = 0.2;//100.0*vMaterialPoint[0]->d_YieldStress/vMaterialPoint[0]->d_ElasticModulus;
//		}
//		if(vMaterialPoint_CPDI.size() != 0)
//		{
//			d6Strain_Minimum[0] = 0.0;//vMaterialPoint_CPDI[0]->d_YieldStress/vMaterialPoint_CPDI[0]->d_ElasticModulus;
//			d6Strain_Maximum[0] = 0.1;//50.0*vMaterialPoint_CPDI[0]->d_YieldStress/vMaterialPoint_CPDI[0]->d_ElasticModulus;
//		}
//		CR.calculateState_J2(d6Strain_Minimum);
//		float fJ2_Minimum = CR.d_J2;
//		CR.calculateState_J2(d6Strain_Maximum);
//		float fJ2_Maximum = CR.d_J2;

		float fJ2_Minimum = 0.0;
		float fJ2_Maximum = 0.1;

		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			glm::vec3 f3Principal = glm::abs(CR.getPrincipal(thisMP->d6_Strain));
			float fJ2 = glm::max(glm::max(f3Principal.x,f3Principal.y),f3Principal.z);
			glm::vec4 f4objectColor = _BLUE;
			if(thisMP->b_Monitor == true)
			{
				if(fJ2 > fJ2_Minimum)
					f4objectColor = (1.0f-fJ2/fJ2_Maximum) * _BLUE + fJ2/fJ2_Maximum * _RED;
				if(fJ2 > fJ2_Maximum)
					f4objectColor  = _GREEN;
			}
			if(thisMP->b3_DisplacementControl.x || thisMP->b3_DisplacementControl.y || thisMP->b3_DisplacementControl.z)
				f4objectColor = _GRAY;

			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		for(int index_MP = 0; index_MP < vMaterialPoint_CPDI.size(); index_MP++)
		{
			MaterialPoint_CPDI_CC *thisMP = vMaterialPoint_CPDI[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			CR.calculateState_J2(thisMP->d6_Strain_Plastic);
			glm::vec3 f3Principal = glm::abs(CR.getPrincipal(thisMP->d6_Strain));
			float fJ2 = glm::max(glm::max(f3Principal.x,f3Principal.y),f3Principal.z);//CR.d_J2;
			glm::vec4 f4objectColor = _BLUE;
			if(fJ2 > fJ2_Minimum)
				f4objectColor = (1.0f-fJ2/fJ2_Maximum) * _BLUE + fJ2/fJ2_Maximum * _GREEN;
			if(fJ2 > fJ2_Maximum)
				f4objectColor  = _RED;
			if(thisMP->b3_DisplacementControl.x || thisMP->b3_DisplacementControl.y || thisMP->b3_DisplacementControl.z)
				f4objectColor = _GRAY;

			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();

			for(int index_Corner = 0; index_Corner < 4; index_Corner++)
			{
				// particle position
				float fSize = 1.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);

				// shadow
				glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
				glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				// camera and model transformation matices
				Transformation glTransformation(thisMP->a_Corner[index_Corner].d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
				glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
				glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
				glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

				gl_Particle_Mesh->Draw();
			}
		}

		// grid points --------------------------------------------------------
		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();
		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(thisGP->b3_Fixed.y == false && thisGP->b3_Fixed.x == false)
				continue;

			// particle position
			float fSize = 5.0*(mpm_PhysicsEngine->d3_Length_World.x) / 1.0e4;
			glm::vec3 f3Size = glm::vec3(fSize,fSize,fSize);
			if(thisGP->b3_Fixed.y == true)
			{
				f3Size.x *= 10.0;
				f3Size.z *= 10.0;
			}
			else if(thisGP->b3_Fixed.x == true)
			{
				f3Size.y *= 10.0;
				f3Size.z *= 10.0;
			}
//				0.002f*glm::vec3(thisGP->b3_Fixed) + glm::vec3(0.00001);
//			Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize,0.01*fSize,fSize));
			Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), f3Size);
			// particle color
			glm::vec4 f4objectColor = _BLACK;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		gl_BasicProgram.unuse();
	}

	if((int)enum_Canvas::J2_STRESS < (int)enum_Canvas::COUNT)
	{// MP parameteric value, von-Mises stress
		v_Canvas_Texture[(int)enum_Canvas::J2_STRESS]->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Diffuse_Texture->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<MaterialPoint_BC *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();
		std::vector<MaterialPoint_CPDI_CC *> vMaterialPoint_CPDI = mpm_PhysicsEngine->getMaterialPoints_CPDI();

		// material points ----------------------------------------------------
		ConstitutiveRelation CR;
		float fJ2_Minimum = 0.0;
		float fJ2_Maximum = 0.0;
		float fJ2_Average = 0.0;
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			fJ2_Average += thisMP->d_Energy_Strain / vMaterialPoint.size();
		}
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			if(thisMP->b_Monitor == false)
				continue;

			float fJ2 = thisMP->d_Energy_Strain;

			//float fJ2 = CR.getState_J2(thisMP->d6_Stress);
//			glm::vec3 f3Principal = glm::abs(CR.getPrincipal(thisMP->d6_Stress));
//			float fJ2 = glm::max(glm::max(f3Principal.x,f3Principal.y),f3Principal.z);

			if(fJ2 > fJ2_Maximum)
				fJ2_Maximum = fJ2;
		}
		fJ2_Maximum = 2.0*fJ2_Average;
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			float fJ2 = thisMP->d_Energy_Strain;
//			float fJ2 = CR.getState_J2(thisMP->d6_Stress);
//			glm::vec3 f3Principal = glm::abs(CR.getPrincipal(thisMP->d6_Stress));
//			float fJ2 = glm::max(glm::max(f3Principal.x,f3Principal.y),f3Principal.z);

			glm::vec4 f4objectColor = f4Color_Default;

			if(thisMP->b3_DisplacementControl.x || thisMP->b3_DisplacementControl.y || thisMP->b3_DisplacementControl.z)
				f4objectColor = f4Color_DisplacementControl;
			else if(thisMP->b_Surface)
				f4objectColor = f4Color_Surface;
			else if(thisMP->b_Monitor == true)
			{
				//float fScale = fJ2/fJ2_Maximum;
				float fScale = (fJ2-fJ2_Minimum)/(fJ2_Maximum-fJ2_Minimum);
				if(fJ2 > fJ2_Minimum)
					f4objectColor = (1.0f-fScale) * f4Color_Default + fScale * f4Color_State_Mid;
				if(fJ2 > fJ2_Maximum)
					f4objectColor  = f4Color_State_Max;
			}
			else if(thisMP->b_Mark_ESO == false)
				f4objectColor = f4Color_ESO_False;

			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		for(int index_MP = 0; index_MP < vMaterialPoint_CPDI.size(); index_MP++)
		{
			MaterialPoint_CPDI_CC *thisMP = vMaterialPoint_CPDI[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			CR.calculateState_J2(thisMP->d6_Stress);
			float fJ2 = CR.d_J2;
//			glm::vec4 f4objectColor = (1.0f-fJ2/fJ2_Maximum) * _BLUE + fJ2/fJ2_Maximum * _RED;

			glm::vec4 f4objectColor = _BLUE;
			if(fJ2 > fJ2_Minimum)
				f4objectColor = (1.0f-fJ2/fJ2_Maximum) * _BLUE + fJ2/fJ2_Maximum * _RED;
			if(fJ2 > fJ2_Maximum)
				f4objectColor  = _RED;

			if(thisMP->b3_DisplacementControl.x || thisMP->b3_DisplacementControl.y || thisMP->b3_DisplacementControl.z)
				f4objectColor = _GRAY;
			if(thisMP->b_Mark_Stress)
				f4objectColor = _GRAY;

			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();

			for(int index_Corner = 0; index_Corner < 4; index_Corner++)
			{
				// particle position
				float fSize = 1.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);

				// shadow
				glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
				glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				// camera and model transformation matices
				Transformation glTransformation(thisMP->a_Corner[index_Corner].d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
				glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
				glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
				glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

				gl_Particle_Mesh->Draw();
			}
		}

		// grid points --------------------------------------------------------
		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();
		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(thisGP->b3_Fixed.y == false && thisGP->b3_Fixed.x == false)
				continue;

			// particle position
			float fSize = 5.0*(mpm_PhysicsEngine->d3_Length_World.x) / 1.0e4;
			glm::vec3 f3Size = glm::vec3(fSize,fSize,fSize);
			if(thisGP->b3_Fixed.y == true)
			{
				f3Size.x *= 10.0;
				f3Size.z *= 10.0;
			}
			else if(thisGP->b3_Fixed.x == true)
			{
				f3Size.y *= 10.0;
				f3Size.z *= 10.0;
			}
//				0.002f*glm::vec3(thisGP->b3_Fixed) + glm::vec3(0.00001);
//			Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize,0.01*fSize,fSize));
			Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), f3Size);
			// particle color
			glm::vec4 f4objectColor = _BLACK;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		gl_BasicProgram.unuse();
	}

	if((int)enum_Canvas::MASSGRADIENT_GP < (int)enum_Canvas::COUNT)
	{// mass gradient, grid points
		v_Canvas_Texture[(int)enum_Canvas::MASSGRADIENT_GP]->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Diffuse_Texture->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		// grid points --------------------------------------------------------
		std::vector<GridPoint *> vGridPoint = mpm_PhysicsEngine->getGridPoints();

		float fMassGradient_Maximum = 0.0;
		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			if(glm::length(thisGP->d3_MassGradient) > fMassGradient_Maximum)
				fMassGradient_Maximum = glm::length(thisGP->d3_MassGradient);
		}

		for(int index_GP = 0; index_GP < vGridPoint.size(); index_GP++)
		{
			GridPoint *thisGP = vGridPoint[index_GP];

			// particle position
			float fSize = 1.0*(mpm_PhysicsEngine->d3_Length_World.x) / 1.0e3;
			glm::vec3 f3Size = glm::vec3(fSize,fSize,fSize);
			// color
			glm::vec4 f4objectColor = _GREEN;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			// nodal positions
//			f4objectColor = _GRAY;
//			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
//			Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), f3Size);
//			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
//			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
//
//			gl_Particle_Mesh->Draw();

			if(glm::length(thisGP->d3_MassGradient) > 0.01*fMassGradient_Maximum)
			{
				// nodal positions
				f4objectColor = _RED;
				if(thisGP->b_Contact_Negative == true && thisGP->b_Contact_Positive == true)
				{
					f4objectColor = _GRAY;
					f3Size *= 5.0f;
				}
				glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
				Transformation glTransformation(thisGP->d3_Position, glm::vec3(0.0, 0.0, 0.0), f3Size);
				glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
				glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

				gl_Particle_Mesh->Draw();

				// nodal gradient vectors
				f4objectColor = _GREEN;
				glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
				Transformation glTransformation_Vector(glm::dvec3(0.0,0.0,0.0), glm::vec3(0.0,0.0,0.0), glm::vec3(1.0));
				m4TransformationMatrix_Model = glTransformation_Vector.GetModelMatrix();
				glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

				glBegin(GL_LINES);
					glm::vec3 f3Start = thisGP->d3_Position;
					glm::vec3 f3End = glm::vec3(thisGP->d3_Position) + 20.0f*fSize/fMassGradient_Maximum * glm::vec3(thisGP->d3_MassGradient);
					glVertex3f(f3Start.x, f3Start.y, f3Start.z);
					glVertex3f(f3End.x, f3End.y, f3End.z);
				glEnd();
			}
		}

		gl_BasicProgram.unuse();
	}

	if((int)enum_Canvas::MASSGRADIENT_MP < (int)enum_Canvas::COUNT)
	{// mass gradient, particles
		v_Canvas_Texture[(int)enum_Canvas::MASSGRADIENT_MP]->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Diffuse_Texture->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<MaterialPoint_BC *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();

		float fMassGradient_Maximum = 0.0;
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			if(glm::length(thisMP->d3_MassGradient) > fMassGradient_Maximum)
				fMassGradient_Maximum = glm::length(thisMP->d3_MassGradient);
		}

		// material points, classic -------------------------------------------
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

//			if(thisMP->b_Surface != true)
//				continue;

			// particle position
			float fSize = 1.0*(mpm_PhysicsEngine->d3_Length_World.x) / 1.0e3;
			glm::vec3 f3Size = glm::vec3(fSize,fSize,fSize);
			// particle color
			glm::vec4 f4objectColor = _BLUE;
			if(thisMP->b_Mark_Stress)
				f4objectColor = _BLUE;
			if(thisMP->b3_DisplacementControl.x || thisMP->b3_DisplacementControl.y || thisMP->b3_DisplacementControl.z)
				f4objectColor = _GRAY;
			if(thisMP->b_Surface)
				f4objectColor = _RED;

			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			//gl_Particle_Mesh->Draw();

			if(glm::length(thisMP->d3_MassGradient) > 0.0*fMassGradient_Maximum)
			{
				// nodal positions
				f4objectColor = _RED;
				glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
				Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), f3Size);
				glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
				glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

				gl_Particle_Mesh->Draw();

				// nodal gradient vectors
				f4objectColor = _GREEN;
				glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
				Transformation glTransformation_Vector(glm::dvec3(0.0,0.0,0.0), glm::vec3(0.0,0.0,0.0), glm::vec3(1.0));
				m4TransformationMatrix_Model = glTransformation_Vector.GetModelMatrix();
				glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

				glBegin(GL_LINES);
					glm::vec3 f3Start = thisMP->d3_Position;
					glm::vec3 f3End = glm::vec3(thisMP->d3_Position) + 20.0f*fSize/fMassGradient_Maximum * glm::vec3(thisMP->d3_MassGradient);
					glVertex3f(f3Start.x, f3Start.y, f3Start.z);
					glVertex3f(f3End.x, f3End.y, f3End.z);
				glEnd();
			}
		}

		gl_BasicProgram.unuse();
	}

	if((int)enum_Canvas::ENERGY_STRAIN < (int)enum_Canvas::COUNT)
	{// MP parameteric value, strain energy
		v_Canvas_Texture[(int)enum_Canvas::ENERGY_STRAIN]->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Diffuse_Texture->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<MaterialPoint_BC *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();
		std::vector<MaterialPoint_CPDI_CC *> vMaterialPoint_CPDI = mpm_PhysicsEngine->getMaterialPoints_CPDI();

//		vMaterialPoint.insert(vMaterialPoint.end(), vMaterialPoint_CPDI.begin(), vMaterialPoint_CPDI.end());
		// material points ----------------------------------------------------
		ConstitutiveRelation CR;
		float fValue_Maximum = 0.5*(vMaterialPoint_CPDI[0]->d_YieldStress * vMaterialPoint_CPDI[0]->d_YieldStress / vMaterialPoint_CPDI[0]->d_ElasticModulus * vMaterialPoint_CPDI[0]->d_Volume_Initial);

		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			float fValue = thisMP->d_Energy_Plastic;
			glm::vec4 f4objectColor = (1.0f-fValue/fValue_Maximum) * _BLUE + fValue/fValue_Maximum * _RED;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		for(int index_MP = 0; index_MP < vMaterialPoint_CPDI.size(); index_MP++)
		{
			MaterialPoint_CPDI_CC *thisMP = vMaterialPoint_CPDI[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			float fValue = thisMP->d_Energy_Strain;
			glm::vec4 f4objectColor = (1.0f-fValue/fValue_Maximum) * _BLUE + fValue/fValue_Maximum * _RED;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();

			for(int index_Corner = 0; index_Corner < 4; index_Corner++)
			{
				// particle position
				float fSize = 1.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);

				// shadow
				glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
				glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				// camera and model transformation matices
				Transformation glTransformation(thisMP->a_Corner[index_Corner].d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
				glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
				glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
				glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

				gl_Particle_Mesh->Draw();
			}
		}

		gl_BasicProgram.unuse();
	}

	if((int)enum_Canvas::ENERGY_PLASTIC < (int)enum_Canvas::COUNT)
	{// MP parameteric value, plastic energy
		v_Canvas_Texture[(int)enum_Canvas::ENERGY_PLASTIC]->bindRenderTarget();
		gl_BasicProgram.use();

		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// diffuse
		GLuint gl_UniformLocation_Diffuse = gl_BasicProgram.getUniformLocation("diffuseTexture");
		glUniform1i(gl_UniformLocation_Diffuse, 0); // save unit 0 for this texture
		gl_Diffuse_Texture->bindTextureUnit(0);
		// shadow
		GLuint gl_UniformLocation_Shadow = gl_BasicProgram.getUniformLocation("shadowTexture");
		glUniform1i(gl_UniformLocation_Shadow, 1); // save unit 0 for this texture
		gl_Shadow_Texture->bindTextureUnit(1);

		GLuint objectColorLocation = gl_BasicProgram.getUniformLocation("objectColor");
		GLuint transformationModelLocation = gl_BasicProgram.getUniformLocation("transformationModelMatrix");
		GLuint transformationCameraLocation = gl_BasicProgram.getUniformLocation("transformationCameraMatrix");
		GLuint transformationShadowLocation = gl_BasicProgram.getUniformLocation("transformationShadowMatrix");
		GLuint lightDirectionLocation = gl_BasicProgram.getUniformLocation("lightDirection");
		GLuint lightColorLocation = gl_BasicProgram.getUniformLocation("lightColor");

		glm::vec3 f3LightDirection = gl_Light->getDirection();
		glm::vec4 f4LightColor = gl_Light->f4_Color;
		glUniform3fv(lightDirectionLocation, 1, &f3LightDirection[0]);
		glUniform4fv(lightColorLocation, 1, &f4LightColor[0]);

		std::vector<MaterialPoint_BC *> vMaterialPoint = mpm_PhysicsEngine->getMaterialPoints();
		std::vector<MaterialPoint_CPDI_CC *> vMaterialPoint_CPDI = mpm_PhysicsEngine->getMaterialPoints_CPDI();

//		vMaterialPoint.insert(vMaterialPoint.end(), vMaterialPoint_CPDI.begin(), vMaterialPoint_CPDI.end());
		// material points ----------------------------------------------------
		ConstitutiveRelation CR;
		float fValue_Maximum = 0.5*(vMaterialPoint_CPDI[0]->d_YieldStress * vMaterialPoint_CPDI[0]->d_YieldStress / vMaterialPoint_CPDI[0]->d_ElasticModulus * vMaterialPoint_CPDI[0]->d_Volume_Initial);

		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			float fValue = thisMP->d_Energy_Plastic;
			glm::vec4 f4objectColor = (1.0f-fValue/fValue_Maximum) * _BLUE + fValue/fValue_Maximum * _RED;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();
		}

		for(int index_MP = 0; index_MP < vMaterialPoint_CPDI.size(); index_MP++)
		{
			MaterialPoint_CPDI_CC *thisMP = vMaterialPoint_CPDI[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			float fValue = thisMP->d_Energy_Plastic;
			glm::vec4 f4objectColor = (1.0f-fValue/fValue_Maximum) * _BLUE + fValue/fValue_Maximum * _RED;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);

			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			Transformation glTransformation(thisMP->d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

			gl_Particle_Mesh->Draw();

			for(int index_Corner = 0; index_Corner < 4; index_Corner++)
			{
				// particle position
				float fSize = 1.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);

				// shadow
				glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
				glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				// camera and model transformation matices
				Transformation glTransformation(thisMP->a_Corner[index_Corner].d3_Position, glm::vec3(0.0, 0.0, 0.0), glm::vec3(fSize));
				glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
				glm::mat4 m4TransformationMatrix_Model = glTransformation.GetModelMatrix();
				glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
				glUniformMatrix4fv(transformationModelLocation, 1, GL_FALSE, &m4TransformationMatrix_Model[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition

				gl_Particle_Mesh->Draw();
			}
		}

		gl_BasicProgram.unuse();
	}

	if(true)
	{// bind the screen for final output
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0); // drawing to the window

		// which .gl program to use
		gl_FinalProgram.use();

		// clear
		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		float fScreenRatio = 1.0 / (int)enum_Canvas::COUNT;
		glm::vec2 f2ScreenRatio = glm::vec2((float)i_ScreenHeight/i_ScreenWidth,1.0);

		if(true)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::MAIN]->bindTextureUnit(0);
			// set viewport
//			float x_Location = (float)enum_Canvas::MAIN / (int)enum_Canvas::COUNT;
//			glViewport(x_Location*i_ScreenWidth, 0.0*i_ScreenHeight, fScreenRatio*i_ScreenWidth, 1.0*i_ScreenHeight);
//			glm::vec2 f2PositionRatio = glm::vec2(0.0/6,0.0);
			glm::vec2 f2PositionRatio = glm::vec2((float)(int)enum_Canvas::MAIN/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
//			gl_Canvas_Mesh->Draw();
			v_Canvas_Mesh[(int)enum_Canvas::MAIN]->Draw();
		}
		if((int)enum_Canvas::SOLID < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::SOLID]->bindTextureUnit(0);
			// set viewport
			glm::vec2 f2PositionRatio = glm::vec2((float)(int)enum_Canvas::SOLID/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::SOLID]->Draw();
		}
		if((int)enum_Canvas::J2_PLASTICSTRAIN < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::J2_PLASTICSTRAIN]->bindTextureUnit(0);
			// set viewport
//			float x_Location = (float)enum_Canvas::J2_PLASTICSTRAIN / (int)enum_Canvas::COUNT;
//			glViewport(x_Location*i_ScreenWidth, 0.0*i_ScreenHeight, fScreenRatio*i_ScreenWidth, 1.0*i_ScreenHeight);
			glm::vec2 f2PositionRatio = glm::vec2((float)(int)enum_Canvas::J2_PLASTICSTRAIN/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::J2_PLASTICSTRAIN]->Draw();
		}
		if((int)enum_Canvas::J2_STRESS < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::J2_STRESS]->bindTextureUnit(0);
			// set viewport
//			float x_Location = (float)enum_Canvas::J2_STRESS / (int)enum_Canvas::COUNT;
//			glViewport(x_Location*i_ScreenWidth, 0.0*i_ScreenHeight, fScreenRatio*i_ScreenWidth, 1.0*i_ScreenHeight);
			glm::vec2 f2PositionRatio = glm::vec2((float)(int)enum_Canvas::J2_STRESS/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::J2_STRESS]->Draw();
		}
		if((int)enum_Canvas::MASSGRADIENT_GP < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::MASSGRADIENT_GP]->bindTextureUnit(0);
			// set viewport
			glm::vec2 f2PositionRatio = glm::vec2((float)(int)enum_Canvas::MASSGRADIENT_GP/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::MASSGRADIENT_GP]->Draw();
		}
		if((int)enum_Canvas::MASSGRADIENT_MP < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::MASSGRADIENT_MP]->bindTextureUnit(0);
			// set viewport
			glm::vec2 f2PositionRatio = glm::vec2((float)(int)enum_Canvas::MASSGRADIENT_MP/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::MASSGRADIENT_MP]->Draw();
		}
		if((int)enum_Canvas::ENERGY_STRAIN < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::ENERGY_STRAIN]->bindTextureUnit(0);
			// set viewport
//			float x_Location = (float)enum_Canvas::ENERGY_STRAIN / (int)enum_Canvas::COUNT;
//			glViewport(x_Location*i_ScreenWidth, 0.0*i_ScreenHeight, fScreenRatio*i_ScreenWidth, 1.0*i_ScreenHeight);
			glm::vec2 f2PositionRatio = glm::vec2((float)(int)enum_Canvas::ENERGY_STRAIN/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::ENERGY_STRAIN]->Draw();
		}
		if((int)enum_Canvas::ENERGY_PLASTIC < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::ENERGY_PLASTIC]->bindTextureUnit(0);
			// set viewport
			glm::vec2 f2PositionRatio = glm::vec2((float)(int)enum_Canvas::ENERGY_PLASTIC/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::ENERGY_PLASTIC]->Draw();
		}

		// undind texture
		gl_FinalProgram.unuse();
	}

	SDL_GL_SwapWindow(p_Window);
}
// ----------------------------------------------------------------------------
