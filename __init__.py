# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from dataclasses import dataclass
from math import sqrt
import math
from typing import List, Tuple
import bpy
import mathutils
from bpy.props import (
    BoolProperty,
    CollectionProperty,
    FloatProperty,
    IntProperty,
    StringProperty,
    EnumProperty,
    PointerProperty,
)
from bpy_extras.io_utils import (
    ImportHelper,
    ExportHelper,
)


bl_info = {
    "name" : "Marble Blast Blender Port",
    "author" : "RandomityGuy",
    "description" : "",
    "blender" : (2, 80, 0),
    "version" : (0, 0, 1),
    "location" : "",
    "warning" : "",
    "category" : "Generic"
}

class PlayMB(bpy.types.Operator):
    """Play MB"""

    bl_idname = "play_mb.mbg"
    bl_label = "Play MB"
    bl_options = {"PRESET"}

    def execute(self, context):
        play_mb()
        return {"FINISHED"}

class MarblePhysicsParams(bpy.types.PropertyGroup):
    enable_mb_physics: BoolProperty(name="Marble Physics", description="Whether to make this object a marble.", default=False)
    marble_diagonal: BoolProperty(name="Diagonal Movement", description="Whether to allow diagonal movement.", default=True)
    marble_speedcap: FloatProperty(name="Speed Cap", description="Maximum speed of marble. (0 to disable)", default=0.0, min=0.0, max=1000.0)
    marble_mass: FloatProperty(name = "Mass", description = "Mass of the marble", default = 1.0, min = 0.1, max = 1000.0)
    marble_max_roll_velocity: FloatProperty(name = "Max Roll Velocity", description = "Maximum roll velocity", default = 15.0, min = 0.1, max = 1000.0)
    marble_angular_acceleration: FloatProperty(name = "Angular Acceleration", description = "Angular acceleration", default = 75.0, min = 0.1, max = 1000.0)
    marble_jump_impulse: FloatProperty(name = "Jump Impulse", description = "Jump impulse", default = 7.5, min = 0.1, max = 1000.0)
    marble_kinetic_friction: FloatProperty(name = "Kinetic Friction", description = "Kinetic friction", default = 0.7, min = 0.1, max = 1000.0)
    marble_static_friction: FloatProperty(name = "Static Friction", description = "Static friction", default = 1.1, min = 0.1, max = 1000.0)
    marble_braking_acceleration: FloatProperty(name = "Braking Acceleration", description = "Braking acceleration", default = 30.0, min = 0.1, max = 1000.0)
    marble_gravity: FloatProperty(name = "Gravity", description = "Gravity", default = 20.0, min = 0.1, max = 1000.0)
    marble_air_accel: FloatProperty(name = "Air Acceleration", description = "Air acceleration", default = 5.0, min = 0.1, max = 1000.0)
    marble_max_dot_slide: FloatProperty(name = "Max Dot Slide", description = "Max dot slide", default = 0.5, min = 0.1, max = 1000.0)
    marble_min_bounce_vel: FloatProperty(name = "Min Bounce Velocity", description = "Min bounce velocity", default = 3.0, min = 0.1, max = 1000.0)
    marble_bounce_kinetic_friction: FloatProperty(name = "Bounce Kinetic Friction", description = "Bounce kinetic friction", default = 0.2, min = 0.1, max = 1000.0)
    marble_bounce_restitution: FloatProperty(name = "Bounce Restitution", description = "Bounce restitution", default = 0.5, min = 0.1, max = 1000.0)

class MarbleMaterialPhysicsParams(bpy.types.PropertyGroup):
    material_friction: FloatProperty(name = "Friction", description = "Friction of the material", default = 1.0, min = -1000.0, max = 1000.0)
    material_restitution: FloatProperty(name = "Restitution", description = "Bounciness of the material", default = 1, min = 0.0, max = 1000.0)
    material_force: FloatProperty(name = "Force", description = "The force this material applies on contact", default = 0, min = 0.0, max = 1000.0)

class MarblePhysicsPanel(bpy.types.Panel):
    bl_label = "Marble Physics Parameters"
    bl_idname = "mb_physicsparams"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "object"

    def draw(self, context):
        layout = self.layout
        obj = context

        sublayout = layout.row()
        sublayout.prop(context.object.mb_props, "enable_mb_physics")
        if context.object.mb_props.enable_mb_physics:
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_diagonal")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_speedcap")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_mass")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_max_roll_velocity")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_angular_acceleration")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_jump_impulse")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_kinetic_friction")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_static_friction")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_braking_acceleration")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_gravity")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_air_accel")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_max_dot_slide")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_min_bounce_vel")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_bounce_kinetic_friction")
            sublayout = layout.row()
            sublayout.prop(context.object.mb_props, "marble_bounce_restitution")

class MarbleMaterialPanel(bpy.types.Panel):
    bl_label = "Marble Material Parameters"
    bl_idname = "mb_materialparams"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "material"

    @classmethod
    def poll(cls, context):
        return (context.material is not None)

    def draw(self, context):
        layout = self.layout
        obj = context.material

        sublayout = layout.row()
        sublayout.prop(obj.mb_props, "material_friction")
        sublayout = layout.row()
        sublayout.prop(obj.mb_props, "material_restitution")
        sublayout = layout.row()
        sublayout.prop(obj.mb_props, "material_force")

class MarblePluginProperties(bpy.types.PropertyGroup):
    camera_sensitivity: FloatProperty(name = "Camera Sensitivity", description = "Camera sensitivity", default = 0.4, min = 0.1, max = 1)

class MarblePluginPanel(bpy.types.Panel):
    bl_label = "Marble Blast"
    bl_idname = "mb_main"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"

    def draw(self, context):
        layout = self.layout
        scene = context.scene

        sublayout = layout.row()
        sublayout.operator(PlayMB.bl_idname, text='Play Marble Blast')
        sublayout = layout.row()
        sublayout.operator(ControlMarble.bl_idname, text='Control Marble')
        sublayout = layout.row()
        sublayout.prop(scene.mb_props, "camera_sensitivity")

controlling_marble = False

input_value = {
    'UP': False,
    'DOWN': False,
    'LEFT': False,
    'RIGHT': False,
    'SPACE': False,
    'mouseX': 0,
    'mouseY': 0
}

input_config = {
    'UP': 'W',
    'DOWN': 'S',
    'LEFT': 'A',
    'RIGHT': 'D',
    'SPACE': 'SPACE'
}

def process_input(event):
    for k, v in input_config.items():
        if event.type == v:
            if event.value == 'PRESS':
                input_value[k] = True
            else:
                input_value[k] = False

    input_value['mouseX'] = event.mouse_x - int(bpy.context.window.width / 2) 
    input_value['mouseY'] = event.mouse_y - int(bpy.context.window.height / 2)

    bpy.context.window.cursor_warp(int(bpy.context.window.width / 2), int(bpy.context.window.height / 2))


class ControlMarble(bpy.types.Operator):
    bl_idname = "control_mb.mbg"
    bl_label = "Control Marble"
    bl_description = "Control Marble with keyboard"

    def invoke(self, context, event):
        global controlling_marble, control_marble
        controlling_marble = True
        if bpy.context.active_object == None:
            return self.report({"ERROR"}, 'Select Marble First.')
        if bpy.context.active_object not in scene_marbles:
            return self.report({"ERROR"}, 'Object is not a Marble!')
        control_marble = scene_marbles[bpy.context.active_object]
        context.window_manager.modal_handler_add(self)
        bpy.context.window.cursor_warp(bpy.context.window.width / 2, bpy.context.window.height / 2)
        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        global controlling_marble
        if event.type == 'ESC':
            controlling_marble = False
            return {'FINISHED'}

        process_input(event)

        return {'RUNNING_MODAL'}

def register():
    bpy.utils.register_class(PlayMB)
    bpy.utils.register_class(ControlMarble)
    bpy.utils.register_class(MarblePhysicsParams)
    bpy.utils.register_class(MarblePhysicsPanel)
    bpy.utils.register_class(MarbleMaterialPhysicsParams)
    bpy.utils.register_class(MarbleMaterialPanel)
    bpy.utils.register_class(MarblePluginProperties)
    bpy.utils.register_class(MarblePluginPanel)

    bpy.types.Object.mb_props = PointerProperty(type=MarblePhysicsParams)
    bpy.types.Material.mb_props = PointerProperty(type=MarbleMaterialPhysicsParams)
    bpy.types.Scene.mb_props = PointerProperty(type=MarblePluginProperties)

def unregister():
    bpy.utils.unregister_class(PlayMB)
    bpy.utils.unregister_class(ControlMarble)
    bpy.utils.unregister_class(MarblePhysicsParams)
    bpy.utils.unregister_class(MarblePhysicsPanel)
    bpy.utils.unregister_class(MarbleMaterialPhysicsParams)
    bpy.utils.unregister_class(MarbleMaterialPanel)
    bpy.utils.unregister_class(MarblePluginProperties)
    bpy.utils.unregister_class(MarblePluginPanel)

    del bpy.types.Object.mb_props
    del bpy.types.Material.mb_props
    del bpy.types.Scene.mb_props

control_marble = None
scene_marbles = {}
scene_obj_positions = {}
scene_obj_velocities = {}

def play_mb():
    bpy.context.scene.render.fps = 60
    bpy.ops.screen.animation_play()
    if tick_mb not in bpy.app.handlers.frame_change_pre:
        bpy.app.handlers.frame_change_pre.append(tick_mb)

    scene_marbles.clear()
    scene_obj_positions.clear()
    scene_obj_velocities.clear()

    for obj in bpy.context.scene.objects:
        if obj.mb_props.enable_mb_physics:
            if obj.type == 'MESH':
                scene_marbles[obj] = Marble(obj)

class Contact:
    def __init__(self, normal: mathutils.Vector, contactDistance: float, obj: bpy.types.Object, restitution: float, friction: float, force: float, velocity: mathutils.Vector = mathutils.Vector((0, 0, 0))):
        self.normal = normal
        self.contactDistance = contactDistance
        self.velocity = velocity
        self.restitution = restitution
        self.friction = friction
        self.force = force
        self.collider = obj

class Move:
    def __init__(self, d: mathutils.Vector, jump: bool, powerup: bool):
        self.d = d
        self.jump = jump
        self.powerup = powerup

class Marble:
    def __init__(self, marbleshape: bpy.types.Object) -> None:
        self.shape = marbleshape
        self.position = marbleshape.location
        self.velocity = mathutils.Vector((0, 0, 0))
        self.omega = mathutils.Vector((0, 0, 0))
        self.diagonal = marbleshape.mb_props.marble_diagonal
        self.speedcap = marbleshape.mb_props.marble_speedcap
        self._mass = marbleshape.mb_props.marble_mass
        self._radius = max(marbleshape.dimensions) / 2
        self._maxRollVelocity = marbleshape.mb_props.marble_max_roll_velocity
        self._angularAcceleration = marbleshape.mb_props.marble_angular_acceleration
        self._jumpImpulse = marbleshape.mb_props.marble_jump_impulse
        self._kineticFriction = marbleshape.mb_props.marble_kinetic_friction
        self._staticFriction = marbleshape.mb_props.marble_static_friction
        self._brakingAcceleration = marbleshape.mb_props.marble_braking_acceleration
        self._gravity = marbleshape.mb_props.marble_gravity
        self._airAccel = marbleshape.mb_props.marble_air_accel
        self._maxDotSlide = marbleshape.mb_props.marble_max_dot_slide
        self._minBounceVel = marbleshape.mb_props.marble_min_bounce_vel
        self._bounceKineticFriction = marbleshape.mb_props.marble_bounce_kinetic_friction
        self._bounceRestitution = marbleshape.mb_props.marble_bounce_restitution
        self._bounceYet = False
        self._bounceSpeed = 0.0
        self._bouncePos = None 
        self._bounceNormal = None
        self._slipAmount = 0.0
        self._contactTime = 0.0
        self._totalTime = 0.0
        self.queuedContacts: List[Contact] = []
        self.gravityDir = mathutils.Vector((0, 0, -1))
        self.contacts: List[Contact] = []
        self.view3d = None
        for a in bpy.context.window.screen.areas:
            if a.type == 'VIEW_3D':
                self.view3d = a
                break

        self.r3d: bpy.types.RegionView3D  = self.view3d.spaces[0].region_3d

    def compute_first_platform_intersect(self, dt):
        shortest_t = dt
        for obj in bpy.context.scene.objects.values():
            if obj == self.shape:
                continue # dont collide with self
            if obj.type != "MESH" or obj.mb_props.enable_mb_physics:
                continue # meshes only
            
            obj_inverse: mathutils.Matrix = obj.matrix_world.inverted_safe()
            obj_inverse_33 = obj_inverse.to_3x3().to_4x4()
            local_vel = obj_inverse_33 @ self.velocity
            result, location, norm, idx = obj.ray_cast(obj_inverse @ self.position, local_vel, distance=(local_vel * shortest_t).length)
            if result:
                shortest_t = min(shortest_t, (location - self.position).length / local_vel.length)

        return max(shortest_t, 0.001)

    def find_contacts(self):
        self.contacts = self.queuedContacts.copy()
        self.queuedContacts.clear()
        mbounds = [[self.position.x - self._radius, self.position.y - self._radius, self.position.z - self._radius], [self.position.x + self._radius, self.position.y + self._radius, self.position.z + self._radius]]
        
        # Static objects
        for obj in bpy.context.scene.objects.values():
            if obj == self.shape:
                continue # dont collide with self
            if obj.type != "MESH" or obj.mb_props.enable_mb_physics:
                continue # meshes only
            # bbox_corners = [obj.matrix_world @ mathutils.Vector(corner) for corner in obj.bound_box]
            # objbounds = [[min([a[0] for a in bbox_corners]), min([a[1] for a in bbox_corners]), min([a[2] for a in bbox_corners])], [max([a[0] for a in bbox_corners]), max([a[1] for a in bbox_corners]), max([a[2] for a in bbox_corners])]]

            # if (mbounds[0][0] > objbounds[1][0] or mbounds[0][1] > objbounds[1][1] or mbounds[0][2] > objbounds[1][2] or mbounds[1][0] < objbounds[0][0] or mbounds[1][0] < objbounds[0][1] or mbounds[2][0] < objbounds[0][2]):
            (found, contactpt, normal, face_index) = obj.closest_point_on_mesh(obj.matrix_world.inverted_safe() @ self.position)

            if found:
                mesh: bpy.types.Mesh = obj.data
                
                contactpt = obj.matrix_world @ mathutils.Vector(contactpt)
                distsq = (self.position - contactpt).length_squared
                if (distsq > self._radius * self._radius):
                    continue

                torquemat = mesh.materials[mesh.polygons[face_index].material_index].mb_props

                dist = self.position - contactpt
                dist.normalize()
                c = Contact(dist, sqrt(distsq), obj, torquemat.material_restitution, torquemat.material_friction, torquemat.material_force, scene_obj_velocities[obj])
                self.contacts.append(c)

        # Marble Objects
        for obj in scene_marbles:
            other_marble = scene_marbles[obj]
            if other_marble == self:
                continue
            distsq = (self.position - other_marble.position).length_squared
            if (distsq > (self._radius + other_marble._radius) * (self._radius + other_marble._radius)):
                continue
            dist = self.position - other_marble.position
            dist.normalize()
            c = Contact(dist, sqrt(distsq), other_marble.shape, 1, 1, 0, other_marble.velocity.copy())
            otherc = Contact(-dist, sqrt(distsq), self.shape, 1, 1, 0, self.velocity.copy())
            self.contacts.append(c)
            # other_marble.queuedContacts.append(otherc)

    def get_marble_axis(self):
        camforward = -mathutils.Vector(self.r3d.view_matrix[2][:3])
        motiondir = (camforward - camforward.project(self.gravityDir))
        motiondir.normalize()
        updir = -self.gravityDir
        sidedir = motiondir.cross(updir)
        sidedir.normalize()
        return [sidedir, motiondir, updir]


    def get_external_forces(self, m: Move, dt: float):
        A = self.gravityDir * self._gravity

        if (len(self.contacts) == 0):
            sideDir, motionDir, upDir = self.get_marble_axis()
            airAccel = self._airAccel
            A = A + ((sideDir * m.d.x) + (motionDir * m.d.y)) * airAccel
        else:
            forceObjCount = 0
            contactForce = 0.0
            contactNormal = mathutils.Vector()
            for contact in self.contacts:
                if contact.force != 0:
                    contactNormal += contact.normal
                    contactForce += contact.force
                    forceObjCount += 1
            if forceObjCount != 0:
                contactNormal.normalize()

                a = contactForce / self._mass
                dot = self.velocity.dot(contactNormal)
                if a > dot:
                    if dot > 0:
                        a -= dot

                    A = A + contactNormal * (a / dt)
                    
        return A

    def compute_move_forces(self, m: Move):
        aControl = mathutils.Vector()
        desiredOmega = mathutils.Vector()
        currentGravityDir = self.gravityDir
        R = currentGravityDir * -self._radius
        rollVelocity = self.omega.cross(R)
        sideDir, motionDir, upDir = self.get_marble_axis()
        currentYVelocity = rollVelocity.dot(motionDir)
        currentXVelocity = rollVelocity.dot(sideDir)
        mv = m.d

        if not self.diagonal:
            mv *= 1.538461565971375
            mvlength = mv.length
            if mvlength > 1:
                mv /= mvlength

        desiredYVelocity = self._maxRollVelocity * mv.y
        desiredXVelocity = self._maxRollVelocity * mv.x

        if (desiredYVelocity != 0 or desiredXVelocity != 0):
            if (currentYVelocity > desiredYVelocity and desiredYVelocity > 0):
                desiredYVelocity = currentYVelocity
            elif (currentYVelocity < desiredYVelocity and desiredYVelocity < 0):
                desiredYVelocity = currentYVelocity

            if (currentXVelocity > desiredXVelocity and desiredXVelocity > 0):
                desiredXVelocity = currentXVelocity
            elif (currentXVelocity < desiredXVelocity and desiredXVelocity < 0):
                desiredXVelocity = currentXVelocity


            rsq = R.length_squared
            desiredOmega = (R.cross(motionDir * desiredYVelocity + sideDir * desiredXVelocity)) * (1 / rsq)
            aControl = desiredOmega - self.omega
            aScalar = aControl.length
            if (aScalar > self._angularAcceleration):
                aControl = aControl * (self._angularAcceleration / aScalar)

            return False, aControl, desiredOmega
        return True, aControl, desiredOmega

    def velocity_cancel(self, surfaceSlide, noBounce):
        SurfaceDotThreshold = 0.0001
        looped = False
        itersIn = 0
        done = False
        while(not done):
            done = True
            itersIn += 1
            for i in range(0, len(self.contacts)):
                sVel = self.velocity - self.contacts[i].velocity
                surfaceDot: float = self.contacts[i].normal.dot(sVel)

                if ((not looped and surfaceDot < 0) or surfaceDot < -SurfaceDotThreshold):
                    velLen = self.velocity.length
                    surfaceVel: mathutils.Vector = self.contacts[i].normal * surfaceDot

                    if (not self._bounceYet):
                        self._bounceYet = True;
                        # playBoundSound(-surfaceDot);
                    

                    if (noBounce):
                        self.velocity = self.velocity - surfaceVel
                    elif (self.contacts[i].collider.mb_props.enable_mb_physics): # is a marble
                        otherMarble: Marble = scene_marbles[self.contacts[i].collider]

                        ourMass = self._mass
                        theirMass = otherMarble._mass

                        bounce = max(self._bounceRestitution, otherMarble._bounceRestitution)

                        # resolve other marble
                        dp = (self.velocity * ourMass) - (self.contacts[i].velocity * theirMass)
                        normP = self.contacts[i].normal * dp.dot(self.contacts[i].normal)

                        normP = normP * (1 + bounce)

                        otherMarble.velocity = self.contacts[i].velocity + (normP * (1 / theirMass))
                        self.contacts[i].velocity = otherMarble.velocity

                        # resolve self
                        normP = self.contacts[i].normal * dp.dot(-self.contacts[i].normal)
                        normP = normP * (1 + bounce)

                        self.velocity = self.velocity + (normP * (1 / ourMass))
                    else:
                        if (self.contacts[i].velocity.length == 0 and not surfaceSlide and surfaceDot > -self._maxDotSlide * velLen):
                            self.velocity = self.velocity - surfaceVel
                            self.velocity.normalize()
                            self.velocity = self.velocity * velLen
                            surfaceSlide = True
                        elif(surfaceDot >= -self._minBounceVel):
                            self.velocity = self.velocity - surfaceVel
                        else:
                            restitution = self._bounceRestitution
                            restitution *= self.contacts[i].restitution

                            velocityAdd = surfaceVel * -(1 + restitution)
                            vAtC = sVel + self.omega.cross(self.contacts[i].normal * -self._radius)
                            normalVel: float = -self.contacts[i].normal.dot(sVel)

                            # bounceEmitter(sVel.length() * restitution, contacts[i].normal);

                            vAtC = vAtC - self.contacts[i].normal * self.contacts[i].normal.dot(sVel)

                            vAtCMag = vAtC.length
                            if (vAtCMag != 0):
                                friction = self._bounceKineticFriction * self.contacts[i].friction

                                angVMagnitude = friction * 5 * normalVel / (2 * self._radius)
                                if (vAtCMag / self._radius < angVMagnitude):
                                    angVMagnitude = vAtCMag / self._radius

                                vAtCDir = vAtC * (1 / vAtCMag)

                                deltaOmega = self.contacts[i].normal.cross(vAtCDir) * angVMagnitude
                                self.omega = self.omega + deltaOmega

                                self.velocity = self.velocity - deltaOmega.cross(self.contacts[i].normal * self._radius)

                            self.velocity += velocityAdd
                    done = False
            looped = True
        if (self.velocity.length_squared < 625):
            gotOne = False;
            dir = mathutils.Vector((0, 0, 0))
            for j in range(0,len(self.contacts)):
                dir2 = dir + self.contacts[j].normal
                if (dir2.length_squared < 0.01):
                    dir2 += self.contacts[j].normal
                dir = dir2
                dir.normalize()
                gotOne = True

            if (gotOne):
                dir.normalize()
                soFar = 0.0
                for k in range(0,len(self.contacts)):
                    dist = self._radius - self.contacts[k].contactDistance;
                    timeToSeparate = 0.1
                    if (dist >= 0):
                        f1 = (self.velocity - self.contacts[k].velocity + dir * soFar).dot(self.contacts[k].normal)
                        f2 = timeToSeparate * f1
                        if (f2 < dist):
                            f3 = (dist - f2) / timeToSeparate
                            soFar += f3 / self.contacts[k].normal.dot(dir)
                if (soFar < -25):
                    soFar = -25
                if (soFar > 25):
                    soFar = 25
                self.velocity = self.velocity + dir * soFar
    
    def apply_contact_forces(self, dt: float, m: Move, isCentered: bool, aControl: mathutils.Vector, desiredOmega: mathutils.Vector, A: mathutils.Vector):
        a = mathutils.Vector((0, 0, 0))
        self._slipAmount = 0
        gWorkGravityDir = self.gravityDir 
        bestSurface = -1
        bestNormalForce = 0.0
        for i in range(0,len(self.contacts)):
            if (True):
                normalForce = -self.contacts[i].normal.dot(A)
                if (normalForce > bestNormalForce):
                    bestNormalForce = normalForce
                    bestSurface = i

        bestContact =  self.contacts[bestSurface] if (bestSurface != -1) else None;
        canJump = bestSurface != -1
        if (canJump and m.jump):
            velDifference = self.velocity - bestContact.velocity
            sv = bestContact.normal.dot(velDifference)
            if (sv < 0):
                sv = 0
            if (sv < self._jumpImpulse):
                self.velocity = self.velocity + bestContact.normal * (self._jumpImpulse - sv)

        for j in range(0, len(self.contacts)):
            normalForce2 = -self.contacts[j].normal.dot(A)
            if (normalForce2 > 0 and self.contacts[j].normal.dot(self.velocity - self.contacts[j].velocity) <= 0.0001):
                A = A + self.contacts[j].normal * normalForce2

        if (bestSurface != -1):
            vAtC = self.velocity + self.omega.cross(bestContact.normal * -self._radius) - (bestContact.velocity)
            vAtCMag = vAtC.length
            slipping = False
            aFriction = mathutils.Vector((0, 0, 0))
            AFriction = mathutils.Vector((0, 0, 0))
            if (vAtCMag != 0):
                slipping = True
                friction = 0.0
                friction = self._kineticFriction * bestContact.friction
                angAMagnitude = 5 * friction * bestNormalForce / (2 * self._radius)
                AMagnitude = bestNormalForce * friction
                totalDeltaV = (angAMagnitude * self._radius + AMagnitude) * dt
                if (totalDeltaV > vAtCMag):
                    fraction = vAtCMag / totalDeltaV
                    angAMagnitude *= fraction
                    AMagnitude *= fraction
                    slipping = False
                vAtCDir = vAtC * (1 / vAtCMag)
                aFriction = -bestContact.normal.cross(-vAtCDir) * angAMagnitude
                AFriction = vAtCDir * -AMagnitude
                self._slipAmount = vAtCMag - totalDeltaV
            if (not slipping):
                R = gWorkGravityDir * -(self._radius)
                aadd = R.cross(A) * (1 / R.length_squared)
                if (isCentered):
                    nextOmega = self.omega + (a * dt)
                    aControl = desiredOmega - nextOmega
                    aScalar = aControl.length
                    if (aScalar > self._brakingAcceleration):
                        aControl = aControl * (self._brakingAcceleration / aScalar)

                Aadd = aControl.cross(bestContact.normal * (-self._radius)) * (-1)
                aAtCMag = (aadd.cross(bestContact.normal * (-self._radius)) + Aadd).length;
                friction2 = 0.0
                friction2 = self._staticFriction * bestContact.friction

                if (aAtCMag > friction2 * bestNormalForce):
                    friction2 = 0
                    friction2 = self._kineticFriction * bestContact.friction
                    Aadd = Aadd * (friction2 * bestNormalForce / aAtCMag)
                A = A + Aadd
                a = a + aadd
            A = A + AFriction
            a = a + aFriction

            lastContactNormal = bestContact.normal
        a = a + aControl
        return [A, a]

    def advance_physics(self, dt: float, m: Move):

        time_remaining = dt
        while time_remaining > 0:
            timestep = 0.008
            if timestep > time_remaining:
                timestep = time_remaining

            timestep = self.compute_first_platform_intersect(timestep)

            time_remaining -= timestep

            self.find_contacts()
            isCentered, aControl, desiredOmega = self.compute_move_forces(m)
            self.velocity_cancel(isCentered, False)
            A = self.get_external_forces(m, timestep)
            A, a = self.apply_contact_forces(timestep, m, isCentered, aControl, desiredOmega, A)
            self.velocity = self.velocity + A * timestep
            self.omega = self.omega + a * timestep
            self.velocity_cancel(isCentered, True)

            vellen = self.velocity.length
            if self.speedcap > 0 and vellen > self.speedcap:
                self.velocity = self.velocity * (self.speedcap / vellen)

            trans, rot, scale = self.shape.matrix_world.decompose()

            self.position = self.position + self.velocity * timestep

            trans = self.position

            deltarot = mathutils.Euler((self.omega.x * timestep, self.omega.y * timestep, self.omega.z * timestep)).to_quaternion()
            deltarot = deltarot @ rot

            finalmat = mathutils.Matrix.LocRotScale(trans, deltarot, scale)

            self.shape.matrix_world = finalmat


def tick_mb(x0, x1):
    DT = 1/60
    m = Move(mathutils.Vector((0, 0, 0)), False, False)
    blankm = Move(mathutils.Vector((0, 0, 0)), False, False)
    if input_value['UP']:
        m.d.y = 1
    if input_value['DOWN']:
        m.d.y = -1
    if input_value['LEFT']:
        m.d.x = -1
    if input_value['RIGHT']:
        m.d.x = 1
    if input_value['SPACE']:
        m.jump = True

    for obj in bpy.context.scene.objects.values():
        if obj.mb_props.enable_mb_physics:
            continue # dont collide with self
        if obj.type != "MESH" or obj.mb_props.enable_mb_physics:
            continue # meshes only
        objpos = obj.location
        prevobjpos = objpos.copy() if obj not in scene_obj_positions else scene_obj_positions[obj]
        scene_obj_positions[obj] = objpos
        scene_obj_velocities[obj] = (objpos - prevobjpos) / DT

    for obj in scene_marbles:
        marble = scene_marbles[obj]
        if control_marble == marble:
            marble.advance_physics(DT, m)
        else:
            marble.advance_physics(DT, blankm)

    if controlling_marble:
        bpy.context.scene.cursor.location = control_marble.position

        for region in (r for r in control_marble.view3d.regions if r.type == 'WINDOW'):
            context_override = {'screen': bpy.context.screen, 'area': control_marble.view3d, 'region': region}
            bpy.ops.view3d.view_center_cursor(context_override)

            r3d: bpy.types.RegionView3D = control_marble.r3d

            view_quat = mathutils.Quaternion(r3d.view_rotation)
            view_eul: mathutils.Euler = view_quat.to_euler()

            cam_sens = 1 / 2500 + (1/100 - 1/2500) * bpy.context.scene.mb_props.camera_sensitivity

            view_eul.z -= input_value['mouseX'] * cam_sens
            view_eul.x += input_value['mouseY'] * cam_sens

            r3d.view_rotation = view_eul.to_quaternion()