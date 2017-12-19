#include "GraphicsEngine.h"

// ----------------------------------------------------------------------------
void GraphicsEngine::drawGame(void)
{
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
			float fSize = 0.5*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			glm::vec4 f4objectColor = _RED;
			if(thisMP->b_Mark_Stress)
				f4objectColor = _GREEN;
			if(thisMP->b_Surface)
			{
				f4objectColor = _BLUE;
				fSize *= 1.0;
			}
			if(thisMP->b_DisplacementControl)
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
			if(thisMP->b_DisplacementControl)
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

			if(thisGP->b3_Fixed == glm::bvec3{false, false, false})	continue;

			// particle position
			glm::vec3 f3Size = glm::vec3(0.0002,0.0002,0.0002);
			if(thisGP->b3_Fixed.y == true)
			{
				f3Size.x = 0.0005;
				f3Size.z = 0.0005;
			}
			else if(thisGP->b3_Fixed.x == true)
			{
				f3Size.y = 0.0005;
				f3Size.z = 0.0005;
			}
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

	if(true)
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
			glm::vec4 f4objectColor = _RED;
			if(thisMP->b_Mark_Stress)
				f4objectColor = _GREEN;
			if(thisMP->b_Surface)
				f4objectColor = _BLACK;
			if(thisMP->b_DisplacementControl)
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
			if(thisMP->b_DisplacementControl)
				f4objectColor = _BLUE;

			if(thisMP->b_DisplacementControl)
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
			float fSize = 0.002;
			glm::vec3 f3Size = glm::vec3(0.0002,0.0002,0.0002);
			if(thisGP->b3_Fixed.y == true)
			{
				f3Size.x = 0.0005;
				f3Size.z = 0.0005;
			}
			else if(thisGP->b3_Fixed.x == true)
			{
				f3Size.y = 0.0005;
				f3Size.z = 0.0005;
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

	if(true)
	{// grid mass gradient
		v_Canvas_Texture[(int)enum_Canvas::MASSGRADIENT]->bindRenderTarget();
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
			float fSize = 0.002;// for problems with dimensions of 1m
			glm::vec3 f3Size = glm::vec3(0.0001,0.0001,0.0001);
			// color
			glm::vec4 f4objectColor = _GREEN;
//			if(thisGP->b_Contact_Positive && thisGP->b_Contact_Negative)
//				f4objectColor = _RED;
//			else
//				continue;
			glUniform4fv(objectColorLocation, 1, &f4objectColor[0]);
			// shadow
			glm::mat4 m4LightTransformation = gl_Light->getViewProjection();
			glUniformMatrix4fv(transformationShadowLocation, 1, GL_FALSE, &m4LightTransformation[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition
			// camera and model transformation matices
			glm::mat4 m4TransformationMatrix_Camera = gl_Camera->getViewProjection();
			glUniformMatrix4fv(transformationCameraLocation, 1, GL_FALSE, &m4TransformationMatrix_Camera[0][0]); // 1 for sending only 1 matrix, GL_FALSE because we don;t want transposition


			if(glm::length(thisGP->d3_MassGradient) > 0.0)
			{
				// nodal positions
				f4objectColor = _RED;
				if(thisGP->b_Contact == true)
				{
					f4objectColor = _GREEN;
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
					glm::vec3 f3End = glm::vec3(thisGP->d3_Position) + fSize/fMassGradient_Maximum * glm::vec3(thisGP->d3_MassGradient);
					glVertex3f(f3Start.x, f3Start.y, f3Start.z);
					glVertex3f(f3End.x, f3End.y, f3End.z);
				glEnd();
			}
		}

		gl_BasicProgram.unuse();
	}

	if(true)
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
		double d6Strain_Minimum[6] = {0, 0, 0, 0, 0, 0};
		double d6Strain_Maximum[6] = {0, 0, 0, 0, 0, 0};
		if(vMaterialPoint.size() != 0)
		{
			d6Strain_Minimum[0] = 0.0;//vMaterialPoint[0]->d_YieldStress/vMaterialPoint[0]->d_ElasticModulus;
			d6Strain_Maximum[0] = 0.1;//100.0*vMaterialPoint[0]->d_YieldStress/vMaterialPoint[0]->d_ElasticModulus;
		}
		if(vMaterialPoint_CPDI.size() != 0)
		{
			d6Strain_Minimum[0] = 0.0;//vMaterialPoint_CPDI[0]->d_YieldStress/vMaterialPoint_CPDI[0]->d_ElasticModulus;
			d6Strain_Maximum[0] = 0.1;//50.0*vMaterialPoint_CPDI[0]->d_YieldStress/vMaterialPoint_CPDI[0]->d_ElasticModulus;
		}
		CR.calculateState_J2(d6Strain_Minimum);
		float fJ2_Minimum = 0.0;//CR.d_J2;
		CR.calculateState_J2(d6Strain_Maximum);
		float fJ2_Maximum = 0.12;//CR.d_J2;

		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			CR.calculateState_J2(thisMP->d6_Strain_Plastic);
			float fJ2 = CR.d_J2;
			glm::vec4 f4objectColor = (1.0f-fJ2/fJ2_Maximum) * _BLUE + fJ2/fJ2_Maximum * _RED;
			if(thisMP->b_DisplacementControl)
				f4objectColor = _GRAY;
			if(thisMP->b_Surface)
				f4objectColor = _BLACK;
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

			if(thisMP->b_DisplacementControl)
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
			float fSize = 0.002;
			glm::vec3 f3Size = glm::vec3(0.0002,0.0002,0.0002);
			if(thisGP->b3_Fixed.y == true)
			{
				f3Size.x = 0.0005;
				f3Size.z = 0.0005;
			}
			else if(thisGP->b3_Fixed.x == true)
			{
				f3Size.y = 0.0005;
				f3Size.z = 0.0005;
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

	if(true)
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
		//float fJ2_Maximum = 1.0/3.0*glm::pow(vMaterialPoint_CPDI[0]->d_YieldStress + 0.0*vMaterialPoint_CPDI[0]->d_Hardening_Isotropic_C1,2);
//		double d6Stress_Minimum[6] = {0.8*vMaterialPoint_CPDI[0]->d_YieldStress, 0, 0, 0, 0, 0};
//		double d6Stress_Maximum[6] = {vMaterialPoint_CPDI[0]->d_YieldStress+0.0057735*vMaterialPoint_CPDI[0]->d_Hardening_Isotropic_C1, 0, 0, 0, 0, 0};
		// for wave speed
		double d6Stress_Minimum[6] = {0.0, 0, 0, 0, 0, 0};
		double d6Stress_Maximum[6] = {0.0, 0, 0, 0, 0, 0};
		if(vMaterialPoint.size() != 0)
		{
			d6Stress_Maximum[0] = vMaterialPoint[0]->p_Material->d_YieldStress + 0.0*vMaterialPoint[0]->d_Hardening_Isotropic_C1;
		}
		if(vMaterialPoint_CPDI.size() != 0)
		{
			d6Stress_Maximum[0] = 1.2*vMaterialPoint_CPDI[0]->p_Material->d_YieldStress;
		}
		CR.calculateState_J2(d6Stress_Maximum);
		float fJ2_Maximum = CR.d_J2;
		CR.calculateState_J2(d6Stress_Minimum);
		float fJ2_Minimum = CR.d_J2;
//		float fJ2_Maximum = 1.0e-12;
//		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
//		{
//			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];
//
//			CR.calculateState_J2(thisMP->d6_Stress);
//			float fJ2 = CR.d_J2;
//
//			if(fJ2 > fJ2_Maximum)
//				fJ2_Maximum = fJ2;
//		}
		for(int index_MP = 0; index_MP < vMaterialPoint.size(); index_MP++)
		{
			MaterialPoint_BC *thisMP = vMaterialPoint[index_MP];

			// particle position
			float fSize = 2.0*0.4*glm::pow(thisMP->d_Volume, 1.0/3.0);
			// particle color
			CR.calculateState_J2(thisMP->d6_Stress);
			float fJ2 = CR.d_J2;
			glm::vec4 f4objectColor = (1.0f-fJ2/fJ2_Maximum) * _BLUE + fJ2/fJ2_Maximum * _RED;
			if(thisMP->b_DisplacementControl)
				f4objectColor = _GRAY;
			if(thisMP->b_Surface)
				f4objectColor = _BLACK;
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

			if(thisMP->b_DisplacementControl)
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
			float fSize = 0.002;
			glm::vec3 f3Size = glm::vec3(0.0002,0.0002,0.0002);
			if(thisGP->b3_Fixed.y == true)
			{
				f3Size.x = 0.0005;
				f3Size.z = 0.0005;
			}
			else if(thisGP->b3_Fixed.x == true)
			{
				f3Size.y = 0.0005;
				f3Size.z = 0.0005;
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

	if(false)
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

	if(false)
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
			glm::vec2 f2PositionRatio = glm::vec2(0.0/6,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			gl_Canvas_Mesh->Draw();
		}
		if((int)enum_Canvas::SOLID < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::SOLID]->bindTextureUnit(0);
			// set viewport
			glm::vec2 f2PositionRatio = glm::vec2(1.0/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::SOLID]->Draw();
		}
		if((int)enum_Canvas::MASSGRADIENT < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::MASSGRADIENT]->bindTextureUnit(0);
			// set viewport
			glm::vec2 f2PositionRatio = glm::vec2(2.0/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::MASSGRADIENT]->Draw();
		}
		if((int)enum_Canvas::J2_PLASTICSTRAIN < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::J2_PLASTICSTRAIN]->bindTextureUnit(0);
			// set viewport
//			float x_Location = (float)enum_Canvas::J2_PLASTICSTRAIN / (int)enum_Canvas::COUNT;
//			glViewport(x_Location*i_ScreenWidth, 0.0*i_ScreenHeight, fScreenRatio*i_ScreenWidth, 1.0*i_ScreenHeight);
			glm::vec2 f2PositionRatio = glm::vec2(3.0/(int)enum_Canvas::COUNT,0.0);
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
			glm::vec2 f2PositionRatio = glm::vec2(4.0/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::J2_STRESS]->Draw();
		}
		if((int)enum_Canvas::ENERGY_STRAIN < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::ENERGY_STRAIN]->bindTextureUnit(0);
			// set viewport
//			float x_Location = (float)enum_Canvas::ENERGY_STRAIN / (int)enum_Canvas::COUNT;
//			glViewport(x_Location*i_ScreenWidth, 0.0*i_ScreenHeight, fScreenRatio*i_ScreenWidth, 1.0*i_ScreenHeight);
			glm::vec2 f2PositionRatio = glm::vec2(5.0/(int)enum_Canvas::COUNT,0.0);
			glViewport(f2PositionRatio.x*i_ScreenWidth, f2PositionRatio.y*i_ScreenHeight, f2ScreenRatio.x*i_ScreenWidth, f2ScreenRatio.y*i_ScreenHeight);
			//draw
			v_Canvas_Mesh[(int)enum_Canvas::ENERGY_STRAIN]->Draw();
		}
		if((int)enum_Canvas::ENERGY_PLASTIC < (int)enum_Canvas::COUNT)
		{
			// bind texture
			v_Canvas_Texture[(int)enum_Canvas::ENERGY_PLASTIC]->bindTextureUnit(0);
			// set viewport
			glm::vec2 f2PositionRatio = glm::vec2(6.0/(int)enum_Canvas::COUNT,0.0);
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
